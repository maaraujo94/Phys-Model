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
  float mass[7]={3.686, 3.556, 3.511, 3.097, 10.355, 10.023, 9.460};
  int stateid = 0;
  
  //read cross section measurements to a vector of form pt/M / cross section*MQuarkonium (nb/GeV) / uncertainty / data ID / K factor / pT/M bin bottom / pT/M bin top / y bin bottom / y bin top / sqrt(s)
  
  //psiprime CMS at 7 TeV
  ns[stateid] = readout(datasigma, "data/CMS_psip_csection.txt", mass[0], 7.9e-3*1000, stateid, 7000, 0);
  stateid++;

  //psiprime ATLAS at 7 TeV
  ns[stateid] = readout(datasigma, "data/ATLAS_psip_csection.txt", mass[0], 5.961*34.46e-4*1000, stateid, 7000, 0);
  stateid++; 
  
  //the ATLAS chic data must have the systematic uncertainty quadratically corrected for the luminosity uncertainty (included as a nuisance parameter)
  //chic2 ATLAS at 7 TeV
  ns[stateid] = readout(datasigma, "data/ATLAS_chic2_csection.txt", mass[1], 1.5*5.961*19.2e-4*1000, stateid, 7000, 1);
  stateid++; 

  //chic1 ATLAS at 7 TeV
  ns[stateid] = readout(datasigma, "data/ATLAS_chic1_csection.txt", mass[2], 1.5*5.961*33.9e-4*1000, stateid, 7000, 1);
  stateid++; 

  //jpsi CMS at 7 TeV
  ns[stateid] = readout(datasigma, "data/CMS_jpsi_csection.txt", mass[3], 5.961e-2*1000, stateid, 7000, 0);
  stateid++; 

  //ups(3S) CMS
  ns[stateid] = readout(datasigma, "data/CMS_ups3_csection.txt", mass[4], 2.4*2.18e4, stateid, 7000, 0);
  stateid++; 

  //ups(3S) ATLAS
  ns[stateid] = readout(datasigma, "data/ATLAS_ups3_csection.txt", mass[4], 2.18e4, stateid, 7000, 0);
  stateid++; 

  //ups(2S) CMS
  ns[stateid] = readout(datasigma, "data/CMS_ups2_csection.txt", mass[5], 2.4*1.93e4, stateid, 7000, 0);
  stateid++; 

  //ups(2S) ATLAS
  ns[stateid] = readout(datasigma, "data/ATLAS_ups2_csection.txt", mass[5], 1.93e4, stateid, 7000, 0);
  stateid++; 

  //ups(1S) CMS
  ns[stateid] = readout(datasigma, "data/CMS_ups1_csection.txt", mass[6], 2.4*2.48e4, stateid, 7000, 0);
  stateid++; 

  //ups(1S) ATLAS
  ns[stateid] = readout(datasigma, "data/ATLAS_ups1_csection.txt", mass[6], 2.48e4, stateid, 7000, 0);
  stateid++; 

  //13 TeV data, has uncertainties in percentage
  //psiprime CMS
  ns[stateid] = readout(datasigma, "data/CMS_psip_csection_13.txt", mass[0], 7.9e-3*1000, stateid, 13000, 2);
  stateid++; 

  //jpsi CMS
  ns[stateid] = readout(datasigma, "data/CMS_jpsi_csection_13.txt", mass[3], 5.961e-2*1000, stateid, 13000, 2);
  stateid++; 

  //ups(3S) CMS
  ns[stateid] = readout(datasigma, "data/CMS_ups3_csection_13.txt", mass[4], 2.18e-2*1000, stateid, 13000, 2);
  stateid++; 

  //ups(2S) CMS
  ns[stateid] = readout(datasigma, "data/CMS_ups2_csection_13.txt", mass[5], 1.93e-2*1000, stateid, 13000, 2);
  stateid++; 

  //ups(1S) CMS
  ns[stateid] = readout(datasigma, "data/CMS_ups1_csection_13.txt", mass[6], 2.48e-2*1000, stateid, 13000, 2);
  stateid++; 

  //cross section ratios
  //chic2/chic1 ratio CMS at 7 TeV
  ns[stateid] = readout(datasigma, "data/CMS_chic_ratio.txt", mass[3], 19.2/33.9, stateid, 7000, 3);
  stateid++; 
  
  //chib2/chib1 ratio CMS
  ns[stateid] = readout(datasigma, "data/CMS_chib_ratio.txt", mass[6], 1., stateid, 7000, 3);
  stateid++; 
  
}
