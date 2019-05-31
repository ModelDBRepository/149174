function Vclamp

figure(1)
clf
grid on

tm = 0:0.1:10;
Concentrations = [0 10e-6 100e-6];
h = zeros(size(Concentrations));
Drug = 'carbamazepine';

protocol = [ ...
	[-inf -50]; ...
	[2 -10]; ...
	[inf -10] ...
	];

junk = zeros(length(tm), length(Concentrations)+1);
junk(:,1) = tm;
for j=1:length(Concentrations)
	concentration = Concentrations(j);
	
	I = INaCK(tm, protocol, concentration, Drug);
	h(1) = line(tm, I, ...
		'linewidth', 2);
	drawnow
	junk(:, j+1) = I;
end

save clampout.dat junk -ascii

str = {};
for j=1:length(Concentrations)
	str{end+1} = sprintf('[%s] = %d{\\mu}M', ...
		Drug, Concentrations(j)*1e6);
end
legend(str)
set(gcf, 'color', 'white')




function I = INaCK(tm, protocol, concentration, drug)
VNa = 50;

q = Q(V(0, protocol), concentration, drug);
[Vx,Dx] = eig(q'); % left eigen vectors are required, hence transpose
[junk indx] = min(abs(diag(Dx))); % find the zero eigen value
P0 = Vx(:,indx); % get the eigen vector
P0 = P0/sum(P0); % Normalise to a probability

[t Y] = ode23(@CKderivs, tm, P0, [], protocol, concentration, drug);
O = Y(:,1);
I = O .* (V(t, protocol)-VNa);

function dPdt = CKderivs(t, P, protocol, concentration, drug)
v = V(t, protocol);
q = Q(v, concentration, drug);
dPdt = (P'*q)';

function I = INaHH(tm, protocol)
VNa = 50;

[minf hinf] = HHrates(protocol.hold, 'inf');

Yinit = [minf hinf];
[t Y] = ode23(@HHderivs, tm, Yinit, [], protocol);

m = Y(:,1);
h = Y(:,2);
I = m.^3 .* h .* (V(t, protocol)-VNa);


function dYdt = HHderivs(t, Y, protocol)
v = V(t, protocol);
m = Y(1);
h = Y(2);

[m_alpha m_beta h_alpha h_beta] = HHrates(v);
dmdt  = m_alpha*(1-m) - m_beta*m;
dhdt  = h_alpha*(1-h) - h_beta*h;

dYdt = [dmdt dhdt]';

function v = V(tm, protocol)
v = zeros(size(tm));
for i=1:length(tm)
	t=tm(i);
	for j=2:size(protocol, 1)
		if protocol(j-1,1)<t && t<=protocol(j,1)
			v(i) = protocol(j-1,2);
			continue
		end
	end
end
