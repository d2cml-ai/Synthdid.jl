include("../src/Synthdid.JL.jl")
california = Synthdid.data("california_prop99")

setup_data = Synthdid.panel_matrices(california)

Y = setup_data.Y;
N0 = setup_data.N0;
T0 = setup_data.T0;

tau_hat = Synthdid.synthdid_estimate(Y, N0, T0);

Synthdid.summary_synth(tau_hat, panel=setup_data);
Synthdid.synthdid_plot(tau_hat)["plot"]