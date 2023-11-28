function [m_a m_b h_a h_b] = HHrates(V, mode)

if nargin==2 && strcmp(mode, 'inf')==1
	m_a = m_inf(V);
	m_b = h_inf(V);
else
	m_a = m_alpha(V);
	m_b = m_beta(V);
	h_a = h_alpha(V);
	h_b = h_beta(V);
end


function x = m_alpha(V)
x = m_inf(V)./tau_m(V);

function x = m_beta(V)
x = (1-m_inf(V))./tau_m(V);

function x = m_inf(V)
Vhalf   = 16.7159;
a       = 10.4440;
x = 1./(1+exp(-(V+Vhalf)/a));

function x = tau_m(V)
tau_m_A = 0.1068;
tau_m_B = 0.0248;
x = tau_m_A*exp(-tau_m_B*V);

function x = h_alpha(V)
x = h_inf(V)./tau_h(V);

function x = h_beta(V)
x = (1-h_inf(V))./tau_h(V);

function x = h_inf(V)
Vhalf_h = 53.6314;
a_h     = -5.5285;
x = 1./(1+exp(-(V+Vhalf_h)/a_h));

function x = tau_h(V)
tau_h_A = 0.4640;
tau_h_B = 0.0712;
x = tau_h_A*exp(-tau_h_B*V);
