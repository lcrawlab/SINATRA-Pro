# SINATRA Pro 

Protein Analysis using Topological Summary Statistics.

see [Github repo](https://github.com/lcrawlab/SINATRA-Pro) for detail.

## Introduction

The sub-image selection problem is to identify physical regions that most explain the variation between two classes of three dimensional shapes. SINATRA is a statistical pipeline for carrying out sub-image analyses using topological summary statistics. The algorithm follows four key steps:

1. 3D shapes of protein structures (represented as triangular meshes) are summarized by a collection of vectors (or curves) detailing their topology (e.g. Euler characteristics, persistence diagrams).
2. A statistical model is used to classify the shapes based on their topological summaries. Here, we make use of a Gaussian process classification model with a probit link function.
3. After fitting the model, an association measure is computed for each topological feature (e.g., centrality measures, posterior inclusion probabilities, p-values, etc).
4. Association measures are mapped back onto the original shapes via a reconstruction algorithm, thus, highlighting evidence of the physical (spatial) locations that best explain the variation between the two groups.

Through detailed simulations, we assess the power of our algorithm as a function of its inputs. As an application of our pipeline, we conduct feature selection for identifying minute conformational changes in five independent protein systems of varying complexities.

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

## Installation and Usage

To install the package, 

        pip3 install SINATRA-Pro

To load the package, 

        import sinatra_pro 

To run the application,

        python3 -m sinatra_pro

## Relevant Citations

Wai Shing Tang*, Gabriel Monteiro da Silva*, Henry Kirveslahti, Erin Skeens, Bibo Feng, Timothy Sudijono, Kevin K. Yang, Sayan Mukherjee, Brenda Rubenstein, and Lorin Crawford. Topological Data Analytic Approach for Discovering Biophysical Signatures in Protein Dynamics.

## License

GNU General Public License v3.0


