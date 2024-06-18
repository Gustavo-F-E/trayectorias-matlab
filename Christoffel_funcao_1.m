function [chr,gg1,gg2,e,g,ev,eu,gv,gu]=Christoffel_funcao_1(S,u,v,ro)
  %Primeira formula fundamental
  e=formula(diff(S,u));
  E1=dot(e,e);
  E=simplify(E1);
  g=formula(diff(S,v));
  G1=dot(g,g);
  G=simplify(G1);
  F1=dot(g,e);
  F=simplify(F1);
  gg1(1,1)=E;
  gg1(1,2)=F;
  gg1(2,1)=F;
  gg1(2,2)=G;
  n=cross(e,g)/norm(cross(e,g));
  L=dot(n,diff(e,u));
  M=dot(n,diff(e,v));
  N=dot(n,diff(g,v));
  gg2(1,1)=L;
  gg2(1,2)=M;
  gg2(2,1)=M;
  gg2(2,2)=N;
  %Símbolos Christoffel
  chr(1,1,1)=simplify((G*(diff(E,u))-2*F*diff(F,u)+F*diff(E,v))/(2*(E*G-F^2)));
  chr(1,1,2)=simplify((2*E*(diff(F,u))-E*diff(E,v)-F*diff(E,u))/(2*(E*G-F^2)));
  chr(1,2,1)=simplify((G*(diff(E,v))-F*diff(G,u))/(2*(E*G-F^2)));
  chr(2,1,1)=chr(1,2,1);
  chr(1,2,2)=simplify((E*(diff(G,u))-F*diff(E,v))/(2*(E*G-F^2)));
  chr(2,1,2)=chr(1,2,2);
  chr(2,2,1)=simplify((2*G*(diff(F,v))-G*diff(G,u)-F*diff(G,v))/(2*(E*G-F^2)));
  chr(2,2,2)=simplify((E*(diff(G,v))-2*F*diff(F,v)-F*diff(G,u))/(2*(E*G-F^2)));
  ev=simplify(diff(gg1(1,1),v));
  eu=simplify(diff(gg1(1,1),u));
  gv=simplify(diff(gg1(2,2),v));
  gu=simplify(diff(gg1(2,2),u));
  end