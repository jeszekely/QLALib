#include <iostream>
#include <cmath>
#include <memory>
#include <algorithm>
#include <vector>

#include "matrix.hpp"
#include "davidson.hpp"

using namespace std;

genMatrix::genMatrix(size_t nr, size_t nc, std::function<vectorMatrix(vectorMatrix&)> h, vector<double> &Diag) :
nrows(nr),
ncols(nc),
H(h)
{
  for (int ii = 0; ii < Diag.size(); ii++)
    diags.push_back(Diag[ii]);//diags(Diag)
}

double genMatrix::diagElem(int ii) {return diags[ii];}

int genMatrix::nc() {return ncols;}

int genMatrix::nr() {return nrows;}

vectorMatrix genMatrix::operator*(vectorMatrix& o){return H(o);}

Davidson::Davidson(genMatrix h, int iVecs, int nVecs, int maxIts, double err) : H(h), initialVecs(iVecs), numVecs(nVecs), maxIterations(maxIts), tolerance(err) {verbose = false;}

tuple<std::shared_ptr<vectorMatrix>,vector<double>> Davidson::diagonalize()
{
  const int n = H.nr(); //Size of a vector

  //eigenvector and eigenvalue return arrays
  auto eigVecs = make_shared<vectorMatrix>(n,numVecs);
  vector<double> eigVals(numVecs, 0.0);

  //initial guess set to unit vectors
  auto trials = make_shared<vectorMatrix>(n, initialVecs);
  for (int ii = 0; ii < initialVecs; ii++)
     trials->element(ii,ii) = H.diagElem(ii);
  trials->makeIdentity();

  //Stored ma
  auto HtStored = make_shared<vectorMatrix>(H**trials);
  int nNewVecs = 0;

  for (int iter = 0; iter < maxIterations; ++iter) {

    //Apply matrix to guess and make subspace matrix, diagonalize
    auto Ht = make_shared<vectorMatrix>(n,HtStored->nc() + nNewVecs);
    Ht->setSub(0,0,*HtStored);
    auto subtrials = make_shared<vectorMatrix>(*trials->getSub(0,HtStored->nc(),n,nNewVecs));
    if (nNewVecs != 0) Ht->setSub(0, HtStored->nc(), H**subtrials);
    HtStored = Ht;

    auto tHt = make_shared<vectorMatrix>(*trials | *Ht);
    auto S = make_shared<vectorMatrix>(*trials | *trials);

    vector<double>S_eigs(S->nc(), 0.0);
    S->diagonalize(S_eigs.data());
    auto tmp = make_shared<vectorMatrix>(S->nc(), count_if(S_eigs.begin(), S_eigs.end(), [] (const double t) { return (t>1.0e-8); }));
    for (int i = 0, current = 0; i < S->nc(); ++i)
      if (S_eigs[i] > 1.0e-8) daxpy_(S->nr(), 1.0/std::sqrt(S_eigs[i]), &S->element(0,i), 1, &tmp->element(0,current++),1);
    vectorMatrix Hprime(*tmp | *tHt * *tmp);

    const int subsize = Hprime.nc();
    vector<double> eigs(subsize, 0.0);
    Hprime.diagonalize(eigs.data());
    vectorMatrix Hsub (*tmp*Hprime);

    vectorMatrix psi( *trials*Hsub );    // current best guesses for eigenvectors
    vectorMatrix sigma( *Ht*Hsub ); // transformed sigma vectors

    copy_n(psi.data(), n*numVecs, eigVecs->data());
    copy_n(eigs.data(), numVecs, eigVals.data());

    vector<shared_ptr<matrixReal>> new_trial_vectors;

    for (int ii = 0; ii < numVecs; ++ii) {
      //matrixReal residual(n, 1);
      //daxpy_(n, -eigs[ii], &psi(0,ii), 1, &residual(0,0), 1);
      //daxpy_(n, 1.0, &sigma(0,ii), 1, &residual(0,0), 1);
      matrixReal residual( (psi.vec(ii)*eigs[ii]) - sigma.vec(ii) );
      if (verbose) cout << residual;

      const double residual_norm = residual.variance();

      if (verbose) cout << setw(6) << iter << setw(6) << ii
                              << fixed << setw(22) << setprecision(12) << eigs[ii]
                              << scientific << setw(16) << setprecision(8) << residual_norm << endl;

      if ( residual_norm > tolerance) {
        auto trial_vector = make_shared<matrixReal>(n, 1);

        // form new guess
        const double en = eigs[ii];
        for (int kk = 0; kk < n; ++kk)
          trial_vector->element(kk, 0) = residual(kk,0) / min(en - H.diagElem(kk), -0.1);

        new_trial_vectors.push_back(trial_vector);
      }
    }
    if (numVecs != 0 && verbose) cout << endl;

    if (new_trial_vectors.empty())
    {
      if (verbose) cout << "Converged!" << endl;
      break;
    }

    auto new_trials = make_shared<vectorMatrix>(n, trials->nc() + new_trial_vectors.size());
    copy_n(trials->data(), trials->size(), new_trials->data());
    for (int ii = 0; ii < new_trial_vectors.size(); ++ii)
      copy_n(new_trial_vectors[ii]->data(), new_trial_vectors[ii]->nr(), &new_trials->element(0, trials->nc() + ii));

    trials = new_trials;
    nNewVecs = new_trial_vectors.size();
  }
  return make_tuple(eigVecs,eigVals);
}
