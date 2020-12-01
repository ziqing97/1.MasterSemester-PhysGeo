% test funktion
transformation=@Meta_transformation;
[L,B] = meshgrid(linspace(-pi,pi,len_l), linspace(-pi/2,pi/2,len_b));
[LM,BM]=transformation(B,L,pi/2 - beta,alpha,gamma);