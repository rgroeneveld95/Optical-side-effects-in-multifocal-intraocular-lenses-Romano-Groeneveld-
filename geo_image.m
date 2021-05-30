function [T1,T2,si2]=geo_image(s01,n_air,n_a,n_iol,n_v,Rc,Ra,Rp,d_iol,V1,V2)

%=============================================================================================
%GEOMETRICAL OPTICS: DETERMINATION OF IMAGE DISTANCE CALCULATION
%=============================================================================================

%SYMBOL EXPLANATION

%Vertices

%V1 = the vertex of the cornea
%V2 = the vertex of the anterior surface of the smooth IOL lens
%V3 = the vertex of the posterior surface of the smooth IOL lens

%Distances 

%So1 = object distance w.r.t. primary principal plane H1
%T1 = distance of the primary principal plane H1 to the primary vertex V1
%d = anterior chamber width measured from vertex V1 to vertex V2
%d_iol = width of the IOL measured from vertex V2 to vertex V3
%T2 = distance of the secondary principal plane H2 to the last vertex V3
%Si2 = image distance w.r.t. secondary principal plane H2

%Radii of curvature

%Rc = corneal curvature
%Ra = anterior surface curvature of the smooth IOL lens
%Rp = posterior surface curvature of the smooth IOL lens

%Refractive indices

%n_air = refractive index of air
%n_a = refractive index aqueous humor
%n_iol = refractive index of the IOL
%n_v = refractive index of the vitreous humor

%=============================================================================================
%CALCULATION
%=============================================================================================

k1=(n_a-n_air)/Rc;
k2=(n_iol-n_a)/Ra;
k3=(n_v-n_iol)/Rp;

Mp=[1 -k2; 0 1];
Md=[1 0;d_iol/n_iol 1];
Ma=[1 -k3; 0 1];

%Transfer matrices of the corneal lens and IOL lens
MV1=[1 -k1;0 1];
MV3V2=Mp*Md*Ma;

%Optical powers of the corneal lens and IOL lens
Pc=-MV1(1,2);
P_iol=-MV3V2(1,2);

%T1 & T2
T1=n_a/MV3V2(1,2)*(MV3V2(1,1)-1);
T2=n_v*(MV3V2(2,2)/MV3V2(1,2)*(MV3V2(1,1)-1)-MV3V2(2,1));

D=V2-V1+T1;

%si2 method 1

% si1=n_a*s01./(Pc*s01-n_air);
% s02=D-si1;
% si2=n_v*s02./(P_iol*s02-n_a);

%si2 method 2

si2=(n_v*(D*Pc-n_a)*s01-n_v*n_air*D)./((P_iol*Pc*D-Pc*n_a-P_iol*n_a)*s01-P_iol*D*n_air+n_a*n_air);
end





