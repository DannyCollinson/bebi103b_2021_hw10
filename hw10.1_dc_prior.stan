functions {
    real lin_growth(real t, real a0, real b) {
        return a0 + b*t;
    }

    real exp_growth(real t, real a0, real k) {
        return a0*e()^(k*t);
    }
}


data {
    // metadata
    int N;
    int J_1;
    int J_2;
    int index_1[J_2];
    int index_2[N];

    // data
    vector[N] time;
    vector[N] area;

    // parameter specifications
    real<lower=0> b_sigma;
    real<lower=0> k_sigma;
    real<lower=0> sigma_b_sigma;
    real<lower=0> sigma_k_sigma;
    real<lower=0> sigma_sigma;
    real<lower=0> a0_mu;
    real<lower=0> a0_sigma;
}


generated quantities {
    vector[N] linear;
    vector[N] exponential;

    real a0_lin = normal_rng(a0_mu, a0_sigma);
    real a0_exp = normal_rng(a0_mu, a0_sigma);
    real b = fabs(normal_rng(0, b_sigma));
    real k = fabs(normal_rng(0, k_sigma));
    real sigma_lin = fabs(normal_rng(0, sigma_sigma));
    real sigma_exp = fabs(normal_rng(0, sigma_sigma));
    real sigma_b = fabs(normal_rng(0, sigma_b_sigma));
    real sigma_k = fabs(normal_rng(0, sigma_k_sigma));

    vector[J_1] b_m;
    vector[J_1] k_m;

    for (i in 1:J_1) {
        b_m[i] = normal_rng(b, sigma_b);
        k_m[i] = normal_rng(k, sigma_k);
    }

    vector[J_2] b_mn;
    vector[J_2] k_mn;

    for (i in 1:J_2) {
        b_mn[i] = normal_rng(b_m[index_1[i]], sigma_b);
        k_mn[i] = normal_rng(k_m[index_1[i]], sigma_k);
    }

    for(i in 1:N) {
        linear[i] = normal_rng(lin_growth(time[i], a0_lin, b_mn[index_2[i]]), sigma_lin);
        exponential[i] = normal_rng(exp_growth(time[i], a0_exp, k_mn[index_2[i]]), sigma_exp);
    }
}