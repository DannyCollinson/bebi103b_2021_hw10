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


parameters {
    real a0_lin_;
    real a0_exp_;
    real<lower=0> b_;
    real<lower=0> k_;
    real<lower=0> sigma_lin_;
    real<lower=0> sigma_exp_;
    real<lower=0> sigma_b_;
    real<lower=0> sigma_k_;

    vector[J_1] b_m_;
    vector[J_1] k_m_;

    vector[J_2] b_mn_;
    vector[J_2] k_mn_;
}


transformed parameters {
    real a0_lin = a0_mu + (a0_lin_ * a0_sigma);
    real a0_exp = a0_mu + (a0_exp_ * a0_sigma);
    real<lower=0> b = b_ * b_sigma;
    real<lower=0> k = k_ * k_sigma;
    real<lower=0> sigma_lin = sigma_lin_ * sigma_sigma;
    real<lower=0> sigma_exp = sigma_exp_ * sigma_sigma;
    real<lower=0> sigma_b = sigma_b_ * sigma_b_sigma;
    real<lower=0> sigma_k = sigma_k_ * sigma_k_sigma;

    vector[J_1] b_m = b + b_m_ * sigma_b;
    vector[J_1] k_m = k + k_m_ * sigma_k;

    vector[J_2] b_mn = b_m[index_1] + b_mn_ * sigma_b;
    vector[J_2] k_mn = k_m[index_1] + k_mn_ * sigma_k;
}


model {
    a0_lin_ ~ normal(0, 1);
    a0_exp_ ~ normal(0, 1);
    b_ ~ normal(0, 1);
    k_ ~ normal(0, 1);
    sigma_lin_ ~ normal(0, 1);
    sigma_exp_ ~ normal(0, 1);
    sigma_b_ ~ normal(0, 1);
    sigma_k_ ~ normal(0, 1);

    b_m_ ~ normal(0, 1);
    k_m_ ~ normal(0, 1);

    b_mn_ ~ normal(0, 1);
    k_mn_ ~ normal(0, 1);

    for (i in 1:N) {
        area[i] ~ normal(lin_growth(time[i], a0_lin, b_mn[index_2[i]]), sigma_lin);
        area[i] ~ normal(exp_growth(time[i], a0_exp, k_mn[index_2[i]]), sigma_exp);
    }
}


generated quantities {
    vector[N] linear;
    vector[N] exponential;

    vector[N] log_lik_lin;
    vector[N] log_lik_exp;

    for(i in 1:N) {
        linear[i] = normal_rng(lin_growth(time[i], a0_lin, b_mn[index_2[i]]), sigma_lin);
        exponential[i] = normal_rng(exp_growth(time[i], a0_exp, k_mn[index_2[i]]), sigma_exp);

        log_lik_lin[i] = normal_lpdf(
            area[i] | lin_growth(time[i], a0_lin, b_mn[index_2[i]]), sigma_lin);
        log_lik_exp[i] = normal_lpdf(
            area[i] | exp_growth(time[i], a0_exp, k_mn[index_2[i]]), sigma_exp);
    }
}