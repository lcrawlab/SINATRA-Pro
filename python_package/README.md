# SINATRA Pro 

Protein Analysis using Topological Summary Statistics.

## Introduction

The sub-image selection problem is to identify physical regions that most explain the variation between two classes of three dimensional shapes. SINATRA is a statistical pipeline for carrying out sub-image analyses using topological summary statistics. The algorithm follows four key steps:

1. 3D shapes of protein structures (represented as triangular meshes) are summarized by a collection of vectors (or curves) detailing their topology (e.g. Euler characteristics, persistence diagrams).
2. A statistical model is used to classify the shapes based on their topological summaries. Here, we make use of a Gaussian process classification model with a probit link function.
3. After fitting the model, an association measure is computed for each topological feature (e.g., centrality measures, posterior inclusion probabilities, p-values, etc).
4. Association measures are mapped back onto the original shapes via a reconstruction algorithm, thus, highlighting evidence of the physical (spatial) locations that best explain the variation between the two groups.

Through detailed simulations, we assess the power of our algorithm as a function of its inputs. Lastly, as an application of our pipeline, we conduct feature selection on a dataset consisting of mandibular molars from different genera of New World Monkeys and examine the physical properties of their teeth that summarize their phylogeny.

## Package Details

Code for implementing the SINATRA Pro pipeline was written in Python 3 (version 3.6.9). As part of this procedure:

1. Reading of trajectory files, alignment of protein structures and neighbor search algorithm are done using the MDAnalysis package (Gowers et. al. 2016, Michaud-Agrawal et. al. 2011).
2. Most athematical calculations are performed using NumPy and SciPy.
3. Inference for the Gaussian process classification (GPC) model was done using elliptical slice sampling (Murray, Prescott, and MacKay 2010).
4. Association measures are computed for the Euler characteristic curves. We use the relative centrality criterion (RATE), which is a variable selection measure for nonlinear and nonparametric statistical methods (Crawford et al. 2019; Ish-Horowicz et al. 2019).

## Dependencies

The SINATRA Pro package depends on these following Python 3 packages.

    numpy >= 1.18.0
    scipy >= 1.5.0
    mdanalysis >= 0.20.0
    fast-histogram >= 0.9
    joblib >= 0.16.0

## Python Package Download

To install the package, 

        pip3 install SINATRA-Pro

To load the package, 

        import sinatra_pro 

To run the application,

        python3 -m sinatra_pro

        usage: __main__.py [-h] [-pa PROTA] [-pb PROTB] [-sa STRUCT_FILE_A]
                           [-ta TRAJ_FILE_A] [-sb STRUCT_FILE_B] [-tb TRAJ_FILE_B]
                           [-dir DIRECTORY] [-pl] [-nc N_CORE] [-n N_SAMPLE]
                           [-of OFFSET] [-s SELECTION] [-r RADIUS] [-hs] [-et EC_TYPE]
                           [-c N_CONE] [-d N_DIRECTION_PER_CONE] [-t CAP_RADIUS]
                           [-l N_FILTRATION] [-bw BANDWIDTH] [-sm SAMPLING_METHOD]
                           [-nm N_MCMC] [-ll] [-v] [-no]

        optional arguments:

              -h, --help            show this help message and exit
              -pa PROTA, --protA PROTA
                                    name of protein A for file naming
              -pb PROTB, --protB PROTB
                                    name of protein B for file naming
              -sa STRUCT_FILE_A, --struct_file_A STRUCT_FILE_A
                                    structure file for protein A (.gro)
              -ta TRAJ_FILE_A, --traj_file_A TRAJ_FILE_A
                                    trajectory file for protein A (.xtc)
              -sb STRUCT_FILE_B, --struct_file_B STRUCT_FILE_B
                                    structure file for protein B (.gro)
              -tb TRAJ_FILE_B, --traj_file_B TRAJ_FILE_B
                                    trajectory file for protein B (.xtc)
              -dir DIRECTORY, --directory DIRECTORY
                                    directory for output files
              -pl, --parallel
              -nc N_CORE, --n_core N_CORE
                                    number of core for parallel computing, default: use
                                    all cores
              -n N_SAMPLE, --n_sample N_SAMPLE
                                    number of sample drawn from trajectory, default: 10
              -of OFFSET, --offset OFFSET
                                    starting frame for sample drawn from trajectory,
                                    default: 0
              -s SELECTION, --selection SELECTION
                                    selection for protein, default: all protein
              -r RADIUS, --radius RADIUS
                                    radius for simplicial construction, default: 2.0
              -hs, --hemisphere     distribute directions over hemisphere instead of whole
                                    sphere
              -et EC_TYPE, --ec_type EC_TYPE
                                    type of Euler characteristic measure (DECT/ECT/SECT),
                                    default: DECT
              -c N_CONE, --n_cone N_CONE
                                    number of cone, default: 1
              -d N_DIRECTION_PER_CONE, --n_direction_per_cone N_DIRECTION_PER_CONE
                                    number of direction per cone, default: 1
              -t CAP_RADIUS, --cap_radius CAP_RADIUS
                                    cap radius, default: 0.8
              -l N_FILTRATION, --n_filtration N_FILTRATION
                                    number of filtration step, default: 20
              -bw BANDWIDTH, --bandwidth BANDWIDTH
                                    bandwidth for elliptical slice sampling, default: 0.01
              -sm SAMPLING_METHOD, --sampling_method SAMPLING_METHOD
                                    sampling method, default: ESS
              -nm N_MCMC, --n_mcmc N_MCMC
                                    number of sample from ESS
              -ll, --logistic_likelihood
                                    use logistic likelihood instead of probit likelihood
              -v, --verbose         verbose
              -no, --name_offset    name folder with offset

## Relevant Citations

Wai Shing Tang*, Gabriel Monteiro da Silva*, Henry Kirveslahti, Erin Skeens, Bibo Feng, Timothy Sudijono, Kevin K. Yang, Sayan Mukherjee, Brenda Rubenstein, and Lorin Crawford. Topological Data Analytic Approach for Discovering Biophysical Signatures in Protein Dynamics.

## License

GNU General Public License v3.0


