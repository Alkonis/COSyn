using namespace arma;

class CollData {

public:
 int layers;
 int steps;

 cube dFdt_0;
 cube tau_0;
 cube dwdn;
 cube Nv;
 
 field<cube> rate_eqtn;
 field<cube> rate_eqtn2;
 field<cube> g;

  CollData(int layers, int steps)  {

  cerr << "Creating cubes." << endl;
  cerr << "Layers,steps:" << endl;
  cerr << layers << endl;
  cerr << steps << endl;
  dFdt_0 = zeros<cube>(10,12,layers);
  tau_0  = zeros<cube>(10,12,layers);
  dwdn   = zeros<cube>(10,12,layers);
  Nv     = zeros<cube>(21,layers,steps);

  cerr << "Creating rate equation." << endl;                                                                                           
  rate_eqtn = field<cube>(steps,1);       
  rate_eqtn2 = field<cube>(steps,1);
  g=field<cube>(10,1);                                                                                                                 
  cerr << "Test1." << endl;                                                                                                            
                                                                                                                                       
  for (int i=0; i<10; i++)                                                                                                             
  {                                                                                                                                    
    g.at(i,0)=zeros<cube>(12,layers,steps);                                                                                            
  }                                                                                                                                    
  cerr <<"Made g." << endl;                                                                                                            
  for (int i=0; i<steps; i++)                                                                           
  {                                                                                                                                    
    cout << i << endl;                                                                                                          
    rate_eqtn.at(i,0)=zeros<cube>(21,21,layers);                                                                                    
    rate_eqtn.at(i,0)(span::all,span(20),span::all).fill(1);                                                                           
  }                                                                                                                                    
/*  for (auto i=0; i<steps;i++) {                                                                                                        
    for (auto j=0; i<layers;j++) {                                                                                                     
      for (auto k=0; k<9;k++) {                                                                                                        
        for (auto q=0; q<10; q++)                                                                                                      
        {                                                                                                                              
          rate_eqtn.at(k+11,0)(q,j,k)=data->einA(k,q);  
        }                                                                                                                              
      }                                                                                                                                
    }                                                                                                                                  
  }  */                                                                                                                                  
                                                                                                                                       
/*  for (int k=0; k<steps; k++)                                                                                                        
 *    {                                                                                                                                    
 *        for (int i=0; i<8 i++)                                                                                                             
 *            {                                                                                                                                  
 *                  rate_eqtn.at(i+1,0).subcube(span(i+1),span::all,span(k))=-vib_einA(i);                                                           
 *                        rate_eqtn.at(i+2,0).subcube(span(i+1,span::all,span(k))=vib_einA(i+1);                                                           
 *                            }
 *                              }*/

  for (auto i=0; i<9; i++)
  {
    //rate_eqtn.at(i+1,0)(span(i+1),span::all,span::all) = vib_einA[i];
        //rate_eqtn.at(i+2,0)(span(i+1),span::all,span::all) = vib_einA[i+1];   //vib_EinA is not defined in this scope--check this and move
  }
   /*rate_eqtn.at(1,0)(span(0),span::all,span::all)=data->vib_einA[0];   //vib_EinA will need to be some kind of Arma class for this to work... I will need to test this
   rate_eqtn.at(10,0)(span(10),span::all,span::all)=data->vib_einA[9];  //look up vib_EinA
   for (auto i=0; i<9; i++)
   {
     rate_eqtn.at(i+11,0)(span(i+1),span::all,span::all) = sum(data->einA.row(i));
   }  */
  //check these threshold-sets ^^^^^^^
  cerr << "Done.";
  
  }
};


