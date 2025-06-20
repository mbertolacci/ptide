#include <Rcpp.h>
#include <RcppEigen.h>
#include <RcppParallel.h>

Eigen::MatrixXd composeCovariance(
    const Eigen::MatrixXd& L,
    const Eigen::VectorXd& sigma
) {
    return (
        sigma.asDiagonal()
        * L
        * L.transpose()
        * sigma.asDiagonal()
    ).selfadjointView<Eigen::Lower>();
}

Eigen::VectorXd vectorSqrt(const Eigen::VectorXd& x) {
    return x.array().sqrt().matrix();
}

//' Calculate covariances between times using a PTide fit
//'
//' Calculate the covariances between times using a PTide fit.
//'
//' @param x A vector of times
//' @param x1 The first vector of times for \code{ptide_cross_covariance}
//' @param x2 A second vector of times for \code{ptide_cross_covariance}
//' @param fit A PTide fit created using \code{\link{ptide_fit}}
//' @param include_noise A vector of 0s and 1s indicating whether to include
//' the white noise term for each time point.
//' @return
//' For \code{ptide_covariance}, a 4-D array of cross-covariances with
//' dimensions n x n x d x d, where n is the length of \code{x} and d is the
//' number of data dimensions. For \code{ptide_cross_covariance}, a 4-D array
//' of cross-covariances with dimensions n1 x n2 x d x d, where n1 and n2 are
//' the lengths of \code{x1} and \code{x2}.
//' @export
// [[Rcpp::export(name = "ptide_covariance", rng=false)]]
Rcpp::NumericVector ptideCovariance(
    const Eigen::VectorXd& x,
    const Rcpp::List& fit,
    const Eigen::VectorXi& include_noise
) {
    using Eigen::Map;
    using Eigen::VectorXd;
    using Eigen::MatrixXd;

    if (x.size() != include_noise.size()) {
        Rcpp::stop("x and include_noise must have the same length");
    }

    Rcpp::List covarianceFunctions = fit["covariance_functions"];
    Rcpp::List firstCovarianceFunction = covarianceFunctions[0];
    VectorXd firstSigmaSquared = firstCovarianceFunction["sigma_squared"];

    int n = x.size();
    int d = firstSigmaSquared.size();
    Rcpp::NumericVector output(n * n * d * d);
    output.attr("dim") = Rcpp::IntegerVector({ n, n, d, d });
    output.fill(0);

    MatrixXd Sigma(n, n);
    for (int i = 0; i < n; ++i) {
        Sigma(i, i) = 1;
    }

    for (int i = 0; i < covarianceFunctions.size(); ++i) {
        Rcpp::List covariance_function = covarianceFunctions[i];
        std::string type = covariance_function["type"];

        VectorXd sigma = vectorSqrt(covariance_function["sigma_squared"]);
        MatrixXd L = covariance_function["L"];

        if (type == "exponential") {
            double lengthDilation = covariance_function["length_dilation"];
            double ell = covariance_function["ell"];
            double ellDilated = ell / lengthDilation;

            for (int col = 0; col < n; ++col) {
                for (int row = col + 1; row < n; ++row) {
                    double d = std::abs(x[row] - x[col]) / ellDilated;
                    Sigma(row, col) = std::exp(-d);
                }
            }
        } else if (type == "matern32") {
            double lengthDilation = covariance_function["length_dilation"];
            double ell = covariance_function["ell"];
            double ellDilated = ell / lengthDilation;

            for (int col = 0; col < n; ++col) {
                for (int row = col + 1; row < n; ++row) {
                    double d = std::sqrt(3) * std::abs(x[row] - x[col]) / ellDilated;
                    Sigma(row, col) = (1.0 + d) * std::exp(-d);
                }
            }
        } else if (type == "matern52") {
            double lengthDilation = covariance_function["length_dilation"];
            double ell = covariance_function["ell"];
            double ellDilated = ell / lengthDilation;

            for (int col = 0; col < n; ++col) {
                for (int row = col + 1; row < n; ++row) {
                    double d = std::sqrt(5) * std::abs(x[row] - x[col]) / ellDilated;
                    Sigma(row, col) = (1.0 + d + d * d / 3.0) * std::exp(-d);
                }
            }
        } else if (type == "squared_exponential") {
            double lengthDilation = covariance_function["length_dilation"];
            double ell = covariance_function["ell"];
            double ellDilated = ell / lengthDilation;

            for (int col = 0; col < n; ++col) {
                for (int row = col + 1; row < n; ++row) {
                    double d = (x[row] - x[col]) / ellDilated;
                    Sigma(row, col) = std::exp(-d * d / 2.0);
                }
            }
        } else if (type == "periodic") {
            double period = covariance_function["period"];
            double ell = covariance_function["ell"];
            for (int col = 0; col < n; ++col) {
                for (int row = col + 1; row < n; ++row) {
                    double distance = std::abs(x[row] - x[col]);
                    double temp = std::sin(M_PI * distance / period);
                    Sigma(row, col) = std::exp(
                        - 2 * temp * temp / (ell * ell)
                    );
                }
            }
        } else if (type == "quasi_periodic") {
            double period = covariance_function["period"];
            double lengthDilation = covariance_function["length_dilation"];
            double ellPeriodic = covariance_function["ell_periodic"];
            double ellDecay = covariance_function["ell_decay"];
            double ellDecayDilated = ellDecay / lengthDilation;
            for (int col = 0; col < n; ++col) {
                for (int row = col + 1; row < n; ++row) {
                    double distance = std::abs(x[row] - x[col]);
                    double temp = std::sin(M_PI * distance / period);
                    Sigma(row, col) = std::exp(
                        - 2 * temp * temp / (ellPeriodic * ellPeriodic)
                        - distance * distance / (2 * ellDecayDilated * ellDecayDilated)
                    );
                }
            }
        }

        MatrixXd Tau = composeCovariance(L, sigma);
        if (type == "white_noise") {
            for (int colD = 0; colD < d; ++colD) {
                for (int rowD = colD; rowD < d; ++rowD) {
                    for (int i = 0; i < n; ++i) {
                        output(
                            i
                            + n * i
                            + n * n * rowD
                            + n * n * d * colD
                        ) += include_noise[i] * Tau(rowD, colD);
                    }
                }
            }
        } else {
            for (int colTau = 0; colTau < d; ++colTau) {
                for (int rowTau = colTau; rowTau < d; ++rowTau) {
                    for (int colSigma = 0; colSigma < n; ++colSigma) {
                        for (int rowSigma = colSigma; rowSigma < n; ++rowSigma) {
                            output(
                                rowSigma
                                + n * colSigma
                                + n * n * rowTau
                                + n * n * d * colTau
                            ) += Tau(rowTau, colTau) * Sigma(rowSigma, colSigma);
                        }
                    }
                }
            }
        }
    }

    // Reflect the outcomes
    for (int colD = 0; colD < d; ++colD) {
        for (int rowD = colD; rowD < d; ++rowD) {
            for (int colN = 0; colN < n; ++colN) {
                for (int rowN = colN; rowN < n; ++rowN) {
                    double value = output(
                        rowN
                        + n * colN
                        + n * n * rowD
                        + n * n * d * colD
                    );
                    output(
                        colN
                        + n * rowN
                        + n * n * rowD
                        + n * n * d * colD
                    ) = value;
                    output(
                        colN
                        + n * rowN
                        + n * n * colD
                        + n * n * d * rowD
                    ) = value;
                    output(
                        rowN
                        + n * colN
                        + n * n * colD
                        + n * n * d * rowD
                    ) = value;
                }
            }
        }
    }

    return output;
}

//' @describeIn ptide_covariance PTide cross covariance
//' @export
// [[Rcpp::export(name = "ptide_cross_covariance", rng=false)]]
Rcpp::NumericVector ptideCrossCovariance(
    const Eigen::VectorXd& x1,
    const Eigen::VectorXd& x2,
    const Rcpp::List& fit
) {
    using Eigen::VectorXd;
    using Eigen::MatrixXd;

    Rcpp::List covarianceFunctions = fit["covariance_functions"];
    Rcpp::List firstCovarianceFunction = covarianceFunctions[0];
    VectorXd firstSigmaSquared = firstCovarianceFunction["sigma_squared"];

    int n1 = x1.rows();
    int n2 = x2.rows();
    int d = firstSigmaSquared.size();
    Rcpp::NumericVector output(n1 * n2 * d * d);
    output.attr("dim") = Rcpp::IntegerVector({ n1, n2, d, d });
    output.fill(0);

    MatrixXd Sigma(n1, n2);

    for (int i = 0; i < covarianceFunctions.size(); ++i) {
        Rcpp::List covariance_function = covarianceFunctions[i];
        std::string type = covariance_function["type"];

        if (type == "white_noise") continue;

        VectorXd sigma = vectorSqrt(covariance_function["sigma_squared"]);
        MatrixXd L = covariance_function["L"];

        if (type == "exponential") {
            double lengthDilation = covariance_function["length_dilation"];
            double ell = covariance_function["ell"];
            double ellDilated = ell / lengthDilation;

            for (int col = 0; col < n2; ++col) {
                for (int row = 0; row < n1; ++row) {
                    double distance = std::abs(x1[row] - x2[col]);
                    Sigma(row, col) = std::exp(-distance / ellDilated);
                }
            }
        } else if (type == "matern32") {
            double lengthDilation = covariance_function["length_dilation"];
            double ell = covariance_function["ell"];
            double ellDilated = ell / lengthDilation;

            for (int col = 0; col < n2; ++col) {
                for (int row = 0; row < n1; ++row) {
                    double d = std::sqrt(3) * std::abs(x1[row] - x2[col]) / ellDilated;
                    Sigma(row, col) = (1.0 + d) * std::exp(-d);
                }
            }
        } else if (type == "matern52") {
            double lengthDilation = covariance_function["length_dilation"];
            double ell = covariance_function["ell"];
            double ellDilated = ell / lengthDilation;

            for (int col = 0; col < n2; ++col) {
                for (int row = 0; row < n1; ++row) {
                    double d = std::sqrt(5) * std::abs(x1[row] - x2[col]) / ellDilated;
                    Sigma(row, col) = (1.0 + d + d * d / 3.0) * std::exp(-d);
                }
            }
        } else if (type == "squared_exponential") {
            double lengthDilation = covariance_function["length_dilation"];
            double ell = covariance_function["ell"];
            double ellDilated = ell / lengthDilation;

            for (int col = 0; col < n2; ++col) {
                for (int row = 0; row < n1; ++row) {
                    double distance = x1[row] - x2[col];
                    Sigma(row, col) = std::exp(
                        - distance * distance / (2 * ellDilated * ellDilated)
                    );
                }
            }
        } else if (type == "periodic") {
            double period = covariance_function["period"];
            double ell = covariance_function["ell"];
            for (int col = 0; col < n2; ++col) {
                for (int row = 0; row < n1; ++row) {
                    double distance = std::abs(x1[row] - x2[col]);
                    double temp = std::sin(M_PI * distance / period);
                    Sigma(row, col) = std::exp(
                        - 2 * temp * temp / (ell * ell)
                    );
                }
            }
        } else if (type == "quasi_periodic") {
            double period = covariance_function["period"];
            double lengthDilation = covariance_function["length_dilation"];
            double ellPeriodic = covariance_function["ell_periodic"];
            double ellDecay = covariance_function["ell_decay"];
            double ellDecayDilated = ellDecay / lengthDilation;
            for (int col = 0; col < n2; ++col) {
                for (int row = 0; row < n1; ++row) {
                    double distance = std::abs(x1[row] - x2[col]);
                    double temp = std::sin(M_PI * distance / period);
                    Sigma(row, col) = std::exp(
                        - 2 * temp * temp / (ellPeriodic * ellPeriodic)
                        - distance * distance / (2 * ellDecayDilated * ellDecayDilated)
                    );
                }
            }
        }

        MatrixXd Tau = composeCovariance(L, sigma);
        for (int colTau = 0; colTau < d; ++colTau) {
            for (int rowTau = 0; rowTau < d; ++rowTau) {
                for (int colSigma = 0; colSigma < n2; ++colSigma) {
                    for (int rowSigma = 0; rowSigma < n1; ++rowSigma) {
                        output(
                            rowSigma
                            + n1 * colSigma
                            + n1 * n2 * rowTau
                            + n1 * n2 * d * colTau
                        ) += Tau(rowTau, colTau) * Sigma(rowSigma, colSigma);
                    }
                }
            }
        }
    }

    return output;
}
