////////////////////////////////////////////////////////////
//auxiliary file that defines the dataset input function  //
////////////////////////////////////////////////////////////

#include <fstream>
using namespace std;

//reads data from file "filename" into vector datasigma,
//according to the auxiliary values given
//methods: 0 - standard, 1 - remove lumi unc from total unc, 2 - unc in %
//         3 - ratio measurement
int readout(vector <vector <double> > &datasigma, string filename, double mass, double signorm, int stateid, double sqrts, int method)
{
  ifstream file;
  string data;
  char nums[12][100];
  double aux = 0.;

  int ns = 0;

  file.open(filename);
  if(file.fail()) {
    cout << "file " << filename << " not found" << endl;
    return 0;
  }
  if(file.is_open())
    while(getline(file,data))
      ns+=1;
  file.clear();
  file.seekg(0);
  for(int j=0; j<ns; j++)
    {
      for(int i=0; i<12; i++)
	file >> nums[i];
      if(method == 1)
	aux=atof(nums[8])*atof(nums[8])-1.8e-2*atof(nums[5])*1.8e-2*atof(nums[5]);
      datasigma[0].push_back(atof(nums[2])/mass);
      if(method == 3) datasigma[1].push_back(atof(nums[5])/signorm);
      else datasigma[1].push_back(atof(nums[5])*mass/signorm);

      if(method == 0)
	datasigma[2].push_back(sqrt(atof(nums[6])*atof(nums[6])+atof(nums[8])*atof(nums[8]))*mass/signorm);
      else if(method == 1)
	datasigma[2].push_back(sqrt(atof(nums[6])*atof(nums[6])+aux)*mass/signorm);
      else if(method == 2)
	datasigma[2].push_back(sqrt(atof(nums[6])*atof(nums[6])*atof(nums[5])*atof(nums[5])+atof(nums[8])*atof(nums[8])*atof(nums[5])*atof(nums[5]))*mass/(100*signorm));
      else if(method == 3)
	datasigma[2].push_back(sqrt(atof(nums[6])*atof(nums[6])+atof(nums[8])*atof(nums[8]))/signorm);

      datasigma[3].push_back(stateid);
      if(method == 3) datasigma[4].push_back(atof(nums[10])-1);
      else datasigma[4].push_back(atof(nums[10])/atof(nums[11])-1);
      datasigma[5].push_back(atof(nums[3])/mass);
      datasigma[6].push_back(atof(nums[4])/mass);
      datasigma[7].push_back(atof(nums[0]));
      datasigma[8].push_back(atof(nums[1]));
      datasigma[9].push_back(sqrts);
    
    }
  
  file.close();

  return ns;
}

// reads data from all desired files, according to the desired corrections
void input(vector <vector <double> > &datasigma, int *ns)
{
  int stateid = 0;
  string fname;
  // folder where I'll read the states
  string orig = "data_first_test/";
  
  //read cross section measurements to a vector of form pt/M / cross section*MQuarkonium (nb/GeV) / uncertainty / data ID / K factor / pT/M bin bottom / pT/M bin top / y bin bottom / y bin top / sqrt(s)
  
  //jpsi CMS at 7 TeV
  ns[stateid] = readout(datasigma, orig+"CMS_jpsi_7_cs.txt", mass[0], 5.961e-2*1e3, stateid, 7000, 0);
  stateid++; 

  //jpsi LHCb at 7 TeV (5 y bins)
  for(int i = 0; i < 5; i++) {
    fname = Form("LHCb_jpsi_7_cs_y%d.txt", i+1);
    fname = orig+fname;
    ns[stateid] = readout(datasigma, fname, mass[0], 1, stateid, 7000, 0);
  stateid++; 
  }
  
  //ups(1S) CMS at 7 TeV
  ns[stateid] = readout(datasigma, orig+"CMS_ups1_7_cs.txt", mass[4], 2.4*2.48e-2*1e6, stateid, 7000, 0);
  stateid++; 

  //ups(1S) ATLAS at 7 TeV (2 y bins)
  for(int i = 0; i < 2; i++) {
    fname = Form("ATLAS_ups1_7_cs_y%d.txt", i+1);
    fname = orig+fname;
    ns[stateid] = readout(datasigma, fname, mass[4], 2.48e-2*1e6, stateid, 7000, 0);
    stateid++; 
  }

  //ups(1S) LHCb at 7 TeV (5 y bins)
  for(int i = 0; i < 5; i++) {
    fname = Form("LHCb_ups1_7_cs_y%d.txt", i+1);
    fname = orig+fname;
    ns[stateid] = readout(datasigma, fname, mass[4], 2.48e-2*1e3, stateid, 7000, 0);
  stateid++; 
  }
    
}
