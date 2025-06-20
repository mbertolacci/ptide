#ifndef STAN_META_HEADER_HPP
#define STAN_META_HEADER_HPP

#include <reduce_sum_vec.hpp>
#include <reduce_sum_vec_prim.hpp>
#include <reduce_sum_vec_dynamic_rev.hpp>
#include <reduce_sum_vec_static_rev.hpp>

struct vecchia_reduce_sum_rsvfunctor {
    template <
        typename T4__,
        typename T5__,
        typename T6__,
        typename T7__,
        typename T8__,
        typename T9__,
        typename T10__,
        typename T11__,
        typename T12__,
        typename T13__,
        typename T14__,
        typename T15__,
        typename T16__,
        typename T17__,
        typename T18__,
        typename T19__,
        typename T20__,
        typename T21__,
        typename T24__,
        typename T25__,
        typename T26__
    >
    Eigen::Matrix<
        stan::promote_args_t<
            T4__, T5__, T6__, T7__, T8__,
            stan::promote_args_t<T9__,
                T10__, T11__, T12__, T13__,
                stan::promote_args_t<
                    T14__, T15__, T16__, T17__, T18__,
                    stan::promote_args_t<
                        T19__, T20__, T21__,
                        stan::value_type_t<T24__>,
                        stan::value_type_t<T25__>,
                        stan::promote_args_t<stan::value_type_t<T26__>>
                    >
                >
            >
        >,
        -1,
        1
    >
    operator()(
        const int& start,
        const int& end,
        std::ostream* pstream__,
        const std::vector<Eigen::Matrix<T4__, -1, 1>>& exponential_sigma,
        const std::vector<Eigen::Matrix<T5__, -1, -1>>& exponential_L,
        const std::vector<T6__>& exponential_ell,
        const std::vector<Eigen::Matrix<T7__, -1, 1>>& matern32_sigma,
        const std::vector<Eigen::Matrix<T8__, -1, -1>>& matern32_L,
        const std::vector<T9__>& matern32_ell,
        const std::vector<Eigen::Matrix<T10__, -1, 1>>& matern52_sigma,
        const std::vector<Eigen::Matrix<T11__, -1, -1>>& matern52_L,
        const std::vector<T12__>& matern52_ell,
        const std::vector<Eigen::Matrix<T13__, -1, 1>>& squared_exponential_sigma,
        const std::vector<Eigen::Matrix<T14__, -1, -1>>& squared_exponential_L,
        const std::vector<T15__>& squared_exponential_ell,
        const std::vector<Eigen::Matrix<T16__, -1, 1>>& quasi_periodic_sigma,
        const std::vector<Eigen::Matrix<T17__, -1, -1>>& quasi_periodic_L,
        const std::vector<T18__>& quasi_periodic_ell_periodic,
        const std::vector<T19__>& quasi_periodic_ell_decay,
        const std::vector<Eigen::Matrix<T20__, -1, 1>>& white_noise_sigma,
        const std::vector<Eigen::Matrix<T21__, -1, -1>>& white_noise_L,
        const std::vector<double>& x,
        const int& d,
        const T24__& y_stacked_arg__,
        const T25__& y_stacked_tilde_arg__,
        const T26__& X_stacked_arg__,
        const std::vector<double>& exponential_length_dilation,
        const std::vector<double>& matern32_length_dilation,
        const std::vector<double>& matern52_length_dilation,
        const std::vector<double>& squared_exponential_length_dilation,
        const std::vector<double>& quasi_periodic_period,
        const std::vector<double>& quasi_periodic_length_dilation,
        const std::vector<int>& block_indices,
        const std::vector<int>& block_last_index,
        const std::vector<int>& block_n_responses,
        const std::vector<int>& block_kernel_group
    ) {
        return vecchia_partial_sums(
            start + 1,
            end,
            exponential_sigma,
            exponential_L,
            exponential_ell,
            matern32_sigma,
            matern32_L,
            matern32_ell,
            matern52_sigma,
            matern52_L,
            matern52_ell,
            squared_exponential_sigma,
            squared_exponential_L,
            squared_exponential_ell,
            quasi_periodic_sigma,
            quasi_periodic_L,
            quasi_periodic_ell_periodic,
            quasi_periodic_ell_decay,
            white_noise_sigma,
            white_noise_L,
            x,
            d,
            y_stacked_arg__,
            y_stacked_tilde_arg__,
            X_stacked_arg__,
            exponential_length_dilation,
            matern32_length_dilation,
            matern52_length_dilation,
            squared_exponential_length_dilation,
            quasi_periodic_period,
            quasi_periodic_length_dilation,
            block_indices,
            block_last_index,
            block_n_responses,
            block_kernel_group,
            pstream__
        );
    }
};

template <
    typename T4__,
    typename T5__,
    typename T6__,
    typename T7__,
    typename T8__,
    typename T9__,
    typename T10__,
    typename T11__,
    typename T12__,
    typename T13__,
    typename T14__,
    typename T15__,
    typename T16__,
    typename T17__,
    typename T18__,
    typename T19__,
    typename T20__,
    typename T21__,
    typename T24__,
    typename T25__,
    typename T26__
>
Eigen::Matrix<
    stan::promote_args_t<
        T4__, T5__, T6__, T7__, T8__,
        stan::promote_args_t<T9__,
            T10__, T11__, T12__, T13__,
            stan::promote_args_t<
                T14__, T15__, T16__, T17__, T18__,
                stan::promote_args_t<
                    T19__, T20__, T21__,
                    stan::value_type_t<T24__>,
                    stan::value_type_t<T25__>,
                    stan::promote_args_t<stan::value_type_t<T26__>>
                >
            >
        >
    >,
    -1,
    1
>
vecchia_reduce_sum(
    const int& start,
    const int& end,
    const int& grain_size,
    const int& strategy,
    const std::vector<Eigen::Matrix<T4__, -1, 1>>& exponential_sigma,
    const std::vector<Eigen::Matrix<T5__, -1, -1>>& exponential_L,
    const std::vector<T6__>& exponential_ell,
    const std::vector<Eigen::Matrix<T7__, -1, 1>>& matern32_sigma,
    const std::vector<Eigen::Matrix<T8__, -1, -1>>& matern32_L,
    const std::vector<T9__>& matern32_ell,
    const std::vector<Eigen::Matrix<T10__, -1, 1>>& matern52_sigma,
    const std::vector<Eigen::Matrix<T11__, -1, -1>>& matern52_L,
    const std::vector<T12__>& matern52_ell,
    const std::vector<Eigen::Matrix<T13__, -1, 1>>& squared_exponential_sigma,
    const std::vector<Eigen::Matrix<T14__, -1, -1>>& squared_exponential_L,
    const std::vector<T15__>& squared_exponential_ell,
    const std::vector<Eigen::Matrix<T16__, -1, 1>>& quasi_periodic_sigma,
    const std::vector<Eigen::Matrix<T17__, -1, -1>>& quasi_periodic_L,
    const std::vector<T18__>& quasi_periodic_ell_periodic,
    const std::vector<T19__>& quasi_periodic_ell_decay,
    const std::vector<Eigen::Matrix<T20__, -1, 1>>& white_noise_sigma,
    const std::vector<Eigen::Matrix<T21__, -1, -1>>& white_noise_L,
    const std::vector<double>& x,
    const int& d,
    const T24__& y_stacked_arg__,
    const T25__& y_stacked_tilde_arg__,
    const T26__& X_stacked_arg__,
    const std::vector<double>& exponential_length_dilation,
    const std::vector<double>& matern32_length_dilation,
    const std::vector<double>& matern52_length_dilation,
    const std::vector<double>& squared_exponential_length_dilation,
    const std::vector<double>& quasi_periodic_period,
    const std::vector<double>& quasi_periodic_length_dilation,
    const std::vector<int>& block_indices,
    const std::vector<int>& block_last_index,
    const std::vector<int>& block_n_responses,
    const std::vector<int>& block_kernel_group,
    std::ostream* pstream__
) {
    int P = X_stacked_arg__.cols();
    if (strategy == 0) {
        return reduce_sum_vec::reduce_sum_vec_dynamic<vecchia_reduce_sum_rsvfunctor>(
            start,
            end,
            2 + P * P + 2 * P,
            grain_size,
            pstream__,
            exponential_sigma,
            exponential_L,
            exponential_ell,
            matern32_sigma,
            matern32_L,
            matern32_ell,
            matern52_sigma,
            matern52_L,
            matern52_ell,
            squared_exponential_sigma,
            squared_exponential_L,
            squared_exponential_ell,
            quasi_periodic_sigma,
            quasi_periodic_L,
            quasi_periodic_ell_periodic,
            quasi_periodic_ell_decay,
            white_noise_sigma,
            white_noise_L,
            x,
            d,
            y_stacked_arg__,
            y_stacked_tilde_arg__,
            X_stacked_arg__,
            exponential_length_dilation,
            matern32_length_dilation,
            matern52_length_dilation,
            squared_exponential_length_dilation,
            quasi_periodic_period,
            quasi_periodic_length_dilation,
            block_indices,
            block_last_index,
            block_n_responses,
            block_kernel_group
        );
    } else {
        return reduce_sum_vec::reduce_sum_vec_static<vecchia_reduce_sum_rsvfunctor>(
            start,
            end,
            2 + P * P + 2 * P,
            grain_size,
            pstream__,
            exponential_sigma,
            exponential_L,
            exponential_ell,
            matern32_sigma,
            matern32_L,
            matern32_ell,
            matern52_sigma,
            matern52_L,
            matern52_ell,
            squared_exponential_sigma,
            squared_exponential_L,
            squared_exponential_ell,
            quasi_periodic_sigma,
            quasi_periodic_L,
            quasi_periodic_ell_periodic,
            quasi_periodic_ell_decay,
            white_noise_sigma,
            white_noise_L,
            x,
            d,
            y_stacked_arg__,
            y_stacked_tilde_arg__,
            X_stacked_arg__,
            exponential_length_dilation,
            matern32_length_dilation,
            matern52_length_dilation,
            squared_exponential_length_dilation,
            quasi_periodic_period,
            quasi_periodic_length_dilation,
            block_indices,
            block_last_index,
            block_n_responses,
            block_kernel_group
        );
    }
}

#endif  // STAN_META_HEADER_HPP
