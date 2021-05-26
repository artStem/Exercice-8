#include <iostream>
#include <fstream>
#include <string>
#include <vector>
#include <valarray>
#include <cmath>
#include <complex> // Pour les nombres complexes
#include "ConfigFile.h"

using namespace std;
typedef vector<complex<double> > vec_cmplx;
long double constexpr M_PI=3.14159265358979323846264338327950288419716939937510582097494459230e0;
// Fonction resolvant le systeme d'equations A * solution = rhs
// ou A est une matrice tridiagonale
template <class T> void triangular_solve(vector<T> const& diag,
                                         vector<T> const& lower,
                                         vector<T> const& upper,
                                         vector<T> const& rhs,
                                         vector<T>& solution)
{
  vector<T> new_diag = diag;
  vector<T> new_rhs = rhs;

  // forward elimination
  for(unsigned int i(1); i<diag.size(); ++i)
  {
    T pivot = lower[i-1] / new_diag[i-1];
    new_diag[i] -= pivot * upper[i-1];
    new_rhs[i] -= pivot * new_rhs[i-1];
  }

  solution.resize(diag.size());

  // solve last equation
  solution[diag.size()-1] = new_rhs[diag.size()-1] / new_diag[diag.size()-1];

  // backward substitution
  for(int i = int(diag.size()) - 2; i >= 0; --i)
  {
    solution[i] = (new_rhs[i] - upper[i] * solution[i+1]) / new_diag[i];
  }
}


// Potentiel V(x) : // TODO
double V(double const& x, double const& omega2, double const& Delta)
{
  return 0.5*omega2*min((x-Delta)*(x-Delta), (x+Delta)*(x+Delta));
}

// Declaration des diagnostics de la particule d'apres sa fonction d'onde psi :
//  - prob calcule la probabilite de trouver la particule entre les points de maillage nL et nR,
//  - E calcule son energie moyenne,
//  - xmoy calcule sa position moyenne,
//  - x2moy calcule sa position au carre moyenne,
//  - pmoy calcule sa quantite de mouvement moyenne,
//  - p2moy calcule sa quantite de mouvement au carre moyenne.
double prob(vec_cmplx const& psi, int nL, unsigned int nR, double dx);
double E(vec_cmplx const& psi, vec_cmplx const& diagH, vec_cmplx const& lowerH, vec_cmplx const& upperH, double const& dx);
double xmoy(vec_cmplx const& psi, const vector<double>& x, double const& dx);
double x2moy(vec_cmplx const& psi, const vector<double>& x, double const& dx);
double pmoy(vec_cmplx const& psi, double const& dx);
double p2moy(vec_cmplx const& psi, double const& dx);

// Fonction pour normaliser une fonction d'onde :
vec_cmplx normalize(vec_cmplx const& psi, double const& dx);

// Les definitions de ces fonctions sont en dessous du main.


int main(int argc,char **argv)
{
  //complex<double> complex_i = complex<double> (0,1); // Nombre imaginaire i

  string inputPath("configuration.in"); // Fichier d'input par defaut
  if(argc>1) // Fichier d'input specifie par l'utilisateur ("./Exercice8 config_perso.in")
    inputPath = argv[1];

  ConfigFile configFile(inputPath); // Les parametres sont lus et stockes dans une "map" de strings.

  for(int i(2); i<argc; ++i) // Input complementaires ("./Exercice8 config_perso.in input_scan=[valeur]")
    configFile.process(argv[i]);

  // Parametres physiques :
  double tfin           = configFile.get<double>("tfin");
  double xL             = configFile.get<double>("xL");
  double xR             = configFile.get<double>("xR");
  double xda             = configFile.get<double>("xda");
  double xdb             = configFile.get<double>("xdb");
  double hbar           = configFile.get<double>("hbar");;
  double m              = configFile.get<double>("mass");
  double omega          = configFile.get<double>("omega");
  double Delta          = configFile.get<double>("Delta");
  double x0             = configFile.get<double>("x0");
  double k0             = 2. * M_PI * double(configFile.get<int>("n")) / (xR-xL);
  double sigma0         = configFile.get<double>("sigma_norm") * (xR-xL);
  double t_detect       = configFile.get<double>("t_detect");
  
  double omega2 = m*omega*omega;

  // Parametres numeriques :
  double dt      = configFile.get<double>("dt");
  int Ninters    = configFile.get<int>("Ninters");
  int Npoints    = Ninters + 1;
  double dx      = (xR-xL) / Ninters;

  // Maillage :
  vector<double> x(Npoints);
  for(int i(0); i<Npoints; ++i)
    x[i] = xL + i*dx;

  // Initialisation de la fonction d'onde :
  vec_cmplx psi(Npoints);
  // TODO: initialiser le paquet d'onde, equation (4.116) du cours
  for(int i(0); i<Npoints; ++i)
    psi[i] =  complex<double> (cos(k0*x[i]),sin(k0*x[i]))*exp(-0.5*pow(((x[i]-x0)/(sigma0)),2.0)); // MODIFY
  // Modifications des valeurs aux bords :
  psi[0] = complex<double> (0.,0.);
  psi[Npoints-1] = complex<double> (0.,0.);
  // Normalisation :
  psi = normalize(psi, dx);

  // Matrices (d: diagonale, a: sous-diagonale, c: sur-diagonale) :
  vec_cmplx dH(Npoints), aH(Ninters), cH(Ninters); // matrice Hamiltonienne
  vec_cmplx dA(Npoints), aA(Ninters), cA(Ninters); // matrice du membre de gauche de l'equation (4.99)
  vec_cmplx dB(Npoints), aB(Ninters), cB(Ninters); // matrice du membre de droite de l'equation (4.99)

  complex<double> a(-hbar*hbar/(2.0*m*dx*dx),0.0); // Coefficient complexe a (devant matrice dérivé seconde)
  complex<double> b(0.0, 0.5*dt/hbar); // Coefficient complexe b (coef multipliant H pour former A et B)

  // TODO: calculer les elements des matrices A, B et H.
  // Ces matrices sont stockees sous forme tridiagonale, d:diagonale, c et a: diagonales superieures et inferieures
  for(int i(0); i<Npoints; ++i) // Boucle sur les points de maillage
  {
    dH[i] = complex<double> (V(x[i], omega2, Delta),0.0) - 2.0*a; // MODIFIER
    dA[i] = 1.0 + b*dH[i]; // MODIFIER
    dB[i] = 1.0 - b*dH[i]; // MODIFIER
  }
  for(int i(0); i<Ninters; ++i) // Boucle sur les intervalles
  {
    aH[i] = cH[i] =  a; // MODIFIER
    aA[i] = cA[i] =  b*cH[i]; // MODIFIER
    aB[i] = cB[i] = -b*cH[i]; // MODIFIER
  }

  // Conditions aux limites: psi nulle aux deux bords
  // TODO: Modifier les matrices A et B pour satisfaire les conditions aux limites
  // INSERER ICI
  dA[0] = dA[Npoints-1] = 1.0;
  dB[0] = dB[Npoints-1] = 0.0;
  cA[0] = aA[0] = cA[Ninters -1] = aA[Ninters-1] = 0.0;
  cB[0] = aB[0] = cB[Ninters -1] = aB[Ninters-1] = 0.0;

  // Fichiers de sortie :
  ofstream fichier_potentiel((configFile.get<string>("output_potential")).c_str());
  fichier_potentiel.precision(15);
  for(int i(0); i<Npoints; ++i)
    fichier_potentiel << x[i] << " " << V(x[i], omega2, Delta) << endl;
  fichier_potentiel.close();

  ofstream fichier_psi((configFile.get<string>("output_squared_wave")).c_str());
  fichier_psi.precision(15);

  ofstream fichier_observables((configFile.get<string>("output_observables")).c_str());
  fichier_observables.precision(15);

  // ecrire position en x
  fichier_psi << 0.0 << " ";
  for(int i(0); i<Npoints-1; ++i){
    fichier_psi << x[i] << " ";
  }
  fichier_psi << x[Npoints-1] << endl;

  // Boucle temporelle :
  valarray<double> print_array=valarray<double>(0.e0,Npoints+1);
  double t,window;
  for(t=0.; t+dt/2.<tfin; t+=dt)
  {
    // Detection de la particule
    if(round(t/dt) == round(t_detect/dt))
    {

      for(int i(0); i<abs(Ninters*xL/(xL-xR)); ++i){
        psi[i] = complex<double> (0.,0.);
      }
      for(int i(abs(Ninters*xL/(xL-xR))); i<abs(Ninters*(xda-xL)/(xL-xR)); ++i){
        window = pow(sin(0.5*M_PI*x[i]/xda),2.e0);
        psi[i] = polar(window*abs(psi[i]),arg(psi[i]));
      }
      for(int i(abs(Ninters*(xdb-xL)/(xL-xR))); i<Ninters; ++i){
        window = pow(0.5*M_PI*x[i]/(xR-xdb),2.0);
        psi[i] = polar(window*abs(psi[i]),arg(psi[i]));
      }
      psi = normalize(psi, dx); // normalise psi pour que la proba totale soit 1
    }

    // Ecriture de |psi|^2 :
    print_array[0] = t;
    for(int i(1); i<Npoints+1; ++i){
      print_array[i] = norm(psi[i-1]); // la fonction C++ norm prend le module au carre
    }
    for(int i(0); i<Npoints; ++i){
      fichier_psi << print_array[i] << " ";
    }
    fichier_psi << print_array[Npoints] << endl;

    // Ecriture des observables :
    fichier_observables << t << " "
                        << prob(psi,0,abs(Ninters*xL/(xL-xR)),dx) << " "              // probabilite que la particule soit en x < 0
                        << prob(psi,abs(Ninters*xL/(xL-xR)),Ninters,dx) << " "        // probabilite que la particule soit en x > 0
                        << prob(psi,0,Ninters,dx) << " " 		   	      // probabilite totale
                        << E(psi,dH,aH,cH,dx) << " "                       	      // Energie
                        << xmoy(psi,x,dx) << " "                           	      // Position moyenne
                        << x2moy(psi,x,dx) << " "                          	      // Position^2 moyenne
                        << pmoy(psi,dx) << " "                             	      // Quantite de mouvement moyenne
                        << p2moy(psi,dx) << " "                           	      // (Quantite de mouvement)^2 moyenne
                        << sqrt(x2moy(psi,x,dx)-xmoy(psi,x,dx)*xmoy(psi,x,dx))*\
			   sqrt(p2moy(psi,dx)-pmoy(psi,dx)*pmoy(psi,dx)) << " "       // Heisenberg index
                        << sqrt(x2moy(psi,x,dx)-xmoy(psi,x,dx)*xmoy(psi,x,dx)) << " " // incertitude en x
                        << sqrt(p2moy(psi,dx)-pmoy(psi,dx)*pmoy(psi,dx)) << endl;     // incertitude en p

    // Calcul du membre de droite :
    vec_cmplx psi_tmp(Npoints,0.);

    // Multiplication psi_tmp = B * psi :
    for(int i(0); i<Npoints; ++i)
      psi_tmp[i] = dB[i] * psi[i];
    for(int i(0); i<Ninters; ++i)
    {
      psi_tmp[i] += cB[i] * psi[i+1];
      psi_tmp[i+1] += aB[i] * psi[i];
    }

    // Resolution de A * psi = psi_tmp :
    triangular_solve(dA, aA, cA, psi_tmp, psi);

  } // Fin de la boucle temporelle

    // ecrire |psi|^2
    print_array[0] = t;
    for(int i(0); i<Npoints; ++i){
      print_array[i+1] = norm(psi[i]);
    }
    for(int i(0); i<Npoints; ++i){
      fichier_psi << print_array[i] << " ";
    }
    fichier_psi << print_array[Npoints] << endl;

    // Ecriture des observables :
    fichier_observables << t << " "
                        << prob(psi,0,abs(Ninters*xL/(xL-xR)),dx) << " "              // probabilite que la particule soit en x < 0
                        << prob(psi,abs(Ninters*xL/(xL-xR)),Ninters,dx) << " "        // probabilite que la particule soit en x > 0
                        << prob(psi,0,Ninters,dx) << " " 		   	      // probabilite totale
                        << E(psi,dH,aH,cH,dx) << " "                       	      // Energie
                        << xmoy(psi,x,dx) << " "                           	      // Position moyenne
                        << x2moy(psi,x,dx) << " "                          	      // Position^2 moyenne
                        << pmoy(psi,dx) << " "                             	      // Quantite de mouvement moyenne
                        << p2moy(psi,dx) << " "                           	      // (Quantite de mouvement)^2 moyenne
                        << sqrt(x2moy(psi,x,dx)-xmoy(psi,x,dx)*xmoy(psi,x,dx))*\
			   sqrt(p2moy(psi,dx)-pmoy(psi,dx)*pmoy(psi,dx)) << " "       // Heisenberg index
                        << sqrt(x2moy(psi,x,dx)-xmoy(psi,x,dx)*xmoy(psi,x,dx)) << " " // incertitude en x
                        << sqrt(p2moy(psi,dx)-pmoy(psi,dx)*pmoy(psi,dx)) << endl;     // incertitude en p

  fichier_observables.close();
  fichier_psi.close();

}


double prob(vec_cmplx const& psi, int nL, unsigned int nR, double dx)
{
  //TODO: calculer la probabilite de trouver la particule entre les points nL et nR
  //Utiliser la formule des trapezes pour l'integration
  double resultat(0.);
  for (unsigned int i(nL); i < nR; ++i) resultat += 0.5*dx*(norm(psi[i]) + norm(psi[i+1])); // MODIFIER
  	// attention : la fonction norm calcule le module carré d'un nombre complexe
  return resultat;
}

//TODO: Calculer les valeurs moyennes des observables E, x, p, x^2, p^2
//Utiliser la formule des trapezes pour l'integration
double E(vec_cmplx const& psi, vec_cmplx const& diagH, vec_cmplx const& lowerH, vec_cmplx const& upperH, double const& dx)
{
  vec_cmplx psi_tmp(psi.size()); // vecteur pour stocker H*psi
  double resultat(0.); // initialiser

  // H(psi): produit de la matrice H et du  vecteur psi
  psi_tmp[0] = diagH[0]*psi[0] + upperH[0]*psi[1];
  for(unsigned int i(1); i<diagH.size()-1; ++i) {
  	psi_tmp[i] = lowerH[i-1]*psi[i-1] + diagH[i]*psi[i] + upperH[i]*psi[i+1];
  }
  psi_tmp.back() = psi[psi.size()-2]*lowerH.back() + psi.back()*diagH.back();
 
  // Integrale de psi* H(psi) dx par la méthode des trapèzes
  for (unsigned int i(0); i<psi.size()-2; ++i) {
  	resultat += 0.5*dx*real(conj(psi[i])*psi_tmp[i] + conj(psi[i+1])*psi_tmp[i+1]);
  } 

  return resultat;
}


double xmoy(vec_cmplx const& psi, const vector<double>& x, double const& dx)
{
  double resultat(0.);
  for (unsigned int i(0); i<psi.size()-2; ++i) {
  	resultat += 0.5*dx*real(conj(psi[i])*x[i]*psi[i] + conj(psi[i+1])*x[i+1]*psi[i+1]);
  }
  return resultat;
}

double x2moy(vec_cmplx const& psi, const vector<double>& x, double const& dx)
{
  double resultat(0.);
  for (unsigned int i(0); i<psi.size()-2; ++i) {
  	resultat += 0.5*dx*real(conj(psi[i])*x[i]*x[i]*psi[i] + conj(psi[i+1])*x[i+1]*x[i+1]*psi[i+1]);
  }
  return resultat;
}


double pmoy(vec_cmplx const& psi, double const& dx)
{
  complex<double> complex_i = complex<double> (0,1); // Nombre imaginaire i
  // unsigned int N(psi.size());
  // Utiliser la definition de p = -i hbar d/dx
  // Utiliser les differences finies centrees pour d/dx
  // Utiliser la formule des trapezes pour l'integration sur x
  // Ignorer la contribution du premier et du dernier point de maillage
  complex<double> resultat(0.,0.);
  // Donner les contributions extremales
  resultat += (conj(psi[0])*(psi[1] - psi[0]));
  resultat += (conj(psi.back())*(psi.back() - psi[psi.size()-2]));
  // Boucler pour avoir les autres contributions
  for (unsigned int i(1); i<psi.size()-2; ++i) resultat += (conj(psi[i])*(psi[i+1] - psi[i-1]));
  // Multiplier par le bon facteur
  resultat *= (-complex_i/2.0);

  return real(resultat);
}


double p2moy(vec_cmplx const& psi, double const& dx)
{
  double resultat(0.);
  // Utiliser la definition de p^2 = - hbar^2 d^2/dx2
  // Utiliser les differences finies centrees pour d^2/dx^2
  // Utiliser la formule des trapezes pour l'integration sur x
  // Ignorer la contribution du premier et du dernier point de maillage
  for (unsigned int i(1); i<psi.size()-2; ++i) resultat += real(conj(psi[i])*(psi[i+1] - 2.0*psi[i] + psi[i-1]));
  resultat /= -dx;
  return resultat;
}


vec_cmplx normalize(vec_cmplx const& psi, double const& dx)
{
  vec_cmplx psi_norm(psi.size());
  double norm = sqrt(prob(psi,0,psi.size()-1,dx));
  for(unsigned int i(0); i<psi.size(); ++i)
    psi_norm[i] = psi[i]/norm;
  return psi_norm;
}





