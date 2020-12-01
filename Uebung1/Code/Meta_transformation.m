% Mit dieser Funktion kann man die Metakoordinaten A,B bekommen
% phi0 und lam0 sind die Breite und Laenge der Metanordpol
% phi und lam sind die Breite und Laenge eines andere Punktes, der am
% Metagreenwich liegt.
% lat und long sind die Breite und Laenge der Punkten, die wir
% transformieren. Die Einheit ist Radiant und der Einsaze als Vektor ist
% erlaubt.
% A,B sind die Metalaenge und Metabreite, die wir brauchen.
% omega0 ist Rotationswinkel.
% 28/05/2019
% Ziqing Yu
function[A,B,omega0]=Meta_transformation(lat,long,phi0,lam0,phi,lam)
    delta_lam=long-lam0;
    omega0=atan2(sin(lam-lam0),sin(phi0)*cos(lam-lam0)-cos(phi0)*tan(phi));
    A=atan2(sin(delta_lam),sin(phi0).*cos(delta_lam)-cos(phi0).*tan(lat))-omega0;
    B=asin(sin(lat)*sin(phi0)+cos(lat).*cos(phi0).*cos(delta_lam));
end

