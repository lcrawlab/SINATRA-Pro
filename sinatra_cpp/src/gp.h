#include <iostream>
#include <armadillo>

using namespace std;
using namespace arma;

//const double pi = M_PI;

dmat calc_covariance_matrix(dmat X,float bandwidth=0.01)
{
    bandwidth = 1.0 / (2.0 * pow(bandwidth,2));
    int n = X.n_cols;    
    dmat K(n,n,fill::eye);
    #pragma omp parallel
    for(int i=0;i<n;i++)
        for(int j=i+1;j<n;j++)
        {
            vec diff = X.col(i) - X.col(j);
            double mean_diff = mean(diff % diff);
            K(i,j) = exp(-mean_diff*bandwidth);
            K(j,i) = K(i,j);
        }
    return K;
}

double probit_log_likelihood(vec latent_variables, vec class_labels)
{
    return accu(log(normcdf(latent_variables % class_labels)));
}

dmat elliptical_slice_sampling(dmat K,vec y,int N_mcmc=100000,int burn_in=1000,int seed=-1,bool verbose=true)
{
    dmat samples;
    int n = K.n_rows;
    int N = y.n_elem;
    if (seed == -1)
        arma_rng::set_seed_random();
    else
        arma_rng::set_seed(seed);
    dmat mcmc_samples(burn_in+N_mcmc,N,fill::zeros);
    dvec mean(n,fill::zeros);
    dmat norm_samples = mvnrnd(mean,K,burn_in+N_mcmc).t();
    dvec unif_samples = randu<dvec>(burn_in+N_mcmc);
    dvec theta = randu<dvec>(burn_in+N_mcmc)*(2.0*pi);
    dvec theta_min = theta - (2.0 * pi);
    dvec theta_max = theta + (2.0 * pi);
    for(int i=1;i<burn_in+N_mcmc;i++)
    {
        if (i < burn_in)
            cout << "Burning in...\r";
        else
            cout << "Elliptical slice sampling Step " << (i-burn_in+1) << "...\r";
        dvec f = conv_to<dvec>::from(mcmc_samples.row(i-1));
        double llh_thresh = probit_log_likelihood(f,y) + log(unif_samples(i));
        dvec f_star = f * cos(theta(i)) + norm_samples.row(i).t();
        while(probit_log_likelihood(f_star,y) < llh_thresh)
        {
            if (theta(i) < 0)
                theta_min(i) = theta(i);
            else
                theta_max(i) = theta(i);
            theta(i) = randu<double>() * (theta_max(i)-theta_min(i)) + theta_min(i);
            f_star = f * cos(theta(i)) + norm_samples.row(i).t() * sin(theta(i));
        }
        mcmc_samples.row(i) = f_star.t();
    }
    cout << '\n';
    return mcmc_samples.rows(burn_in,burn_in+N_mcmc-1);
}

dmat sherman_r(dmat A, dvec u, dvec v)
{
    dmat x = v.t() * A * u;
    dmat B = A * u * v.t() * A;
    double c = 1.0 / (x(0) + 1.0);
    return A - B * c;
}

dvec calc_rate(dmat X, dmat f_draws, bool verbose=true)
{
    dmat beta_draws = (pinv(X) * f_draws.t()).t();
    dmat V = cov(beta_draws);
    dmat D = pinv(V);
    dmat D_U, D_V;
    dvec D_s;
    svd(D_U,D_s,D_V,D);
    uvec ind = find(D_s > 1e-10);
    D_U = D_U.cols(ind);
    dmat U = D_U.each_row() % sqrt(D_s.elem(ind)).t();
    dvec mu = conv_to<dvec>::from(mean(beta_draws,0));
    mu = abs(mu);
    dmat Lambda = U * U.t();
    dvec kld(mu.n_elem,fill::zeros);
    for(int q=0;q<mu.n_elem;q++)
    {
        if (verbose)
            cout << "Calculating KLD(" << q << ")...\r";
        dvec Vq = conv_to<dvec>::from(V.col(q));
        dmat U_Lambda_sub = sherman_r(Lambda,Vq,Vq);
        dmat U_no_q = U_Lambda_sub;
        U_no_q.shed_row(q);
        dmat U_no_qq = U_no_q;
        U_no_qq.shed_col(q);
        dvec U_no_q_q = conv_to<dvec>::from(U_no_q.col(q));
        dvec alpha = U_no_q_q.t() * U_no_qq * U_no_q_q;
        kld(q) = pow(mu(q),2.0) * alpha(0) * 0.5;
    }
    kld.save("kld.txt",raw_ascii);

    if (verbose)
        cout << "KLD calculation Completed.\n";
    
    // Compute the corresponding “RelATive cEntrality” (RATE) measure
    dvec rates = kld / sum(kld);

    /// Find the entropic deviation from a uniform distribution 
    //delta = np.sum(rates*np.log(len(mu)*rates))

    // Calibrate Delta via the effective sample size (ESS) measures from importance sampling ###
    // (Gruber and West, 2016, 2017)
    //eff_samp_size = 1./(1.+delta)*100.
    
    return rates;
}

dmat find_rate_variables_with_other_sampling_methods(dmat X,vec y,float bandwidth = 0.01, string sampling_method = "ESS", int size = 100000, int N_mcmc = 100000, int burn_in = 1000, bool probit = true, int seed = -1)
{   
    int n_fils = X.n_cols;
    uvec nonzero_col = find(any(X,0));
    X = conv_to<dmat>::from(X.cols(nonzero_col));
    cout << "X " << X.n_rows << ' ' << X.n_cols << endl;
    dmat X_colmean = mean(X,0);
    dmat X_colstd = stddev(X,0,0);    
    X.each_row() -= X_colmean;
    X.each_row() /= X_colstd;
    
    cout << X.n_rows << ' ' << X.n_cols << endl;
    dmat K = calc_covariance_matrix(X.t(),bandwidth);
    //K.save("K.bin",arma_binary);
    dmat samples;
    //if (not sampling_method.compare("ESS"))
    samples = elliptical_slice_sampling(K,y,N_mcmc,burn_in,seed,true);
    //samples.save("ESS.bin",arma_binary);    
    //samples.load("ESS.bin",arma_binary);
    dvec rates_nv = calc_rate(X,samples);
    dvec rates(n_fils,fill::zeros);
    for(int i=0;i<nonzero_col.n_elem;i++)
        rates(nonzero_col(i)) = rates_nv(i);
    return rates;
}
