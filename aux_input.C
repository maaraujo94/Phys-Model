#include <fstream>
using namespace std;

void input(vector <vector <double> > &datasigma, int *ns)
{
  float mass[7]={3.686, 3.556, 3.511, 3.097, 10.355, 10.023, 9.460};

  ifstream file;
  string data;
  char nums[12][100];
  int stateid = 0;
  double aux;
  
  //read cross section measurements to a vector of form pt/M / cross section*MQuarkonium (nb/GeV) / uncertainty / data ID / K factor / pT/M bin bottom / pT/M bin top / y bin bottom / y bin top / sqrt(s)/M
  
  //psiprime CMS at 7 TeV
  file.open("data/CMS_psip_csection.txt");
  if(file.is_open())
    while(getline(file,data))
      ns[stateid]+=1;
  file.clear();
  file.seekg(0);
  for(int j=0; j<ns[stateid]; j++)
    {
      for(int i=0; i<12; i++)
	file >> nums[i];
      datasigma[0].push_back(atof(nums[2])/mass[0]);
      datasigma[1].push_back(atof(nums[5])*mass[0]/(7.9e-3*1000));
      datasigma[2].push_back(sqrt(atof(nums[6])*atof(nums[6])+atof(nums[8])*atof(nums[8]))*mass[0]/(7.9e-3*1000));
      datasigma[3].push_back(stateid);
      datasigma[4].push_back(atof(nums[10])/atof(nums[11])-1);
      datasigma[5].push_back(atof(nums[3])/mass[0]);
      datasigma[6].push_back(atof(nums[4])/mass[0]);
      datasigma[7].push_back(atof(nums[0]));
      datasigma[8].push_back(atof(nums[1]));
      datasigma[9].push_back(7000);
    }
  file.close();
  stateid++;

  //psiprime ATLAS at 7 TeV
  file.open("data/ATLAS_psip_csection.txt");
  if(file.is_open())
    while(getline(file,data))
      ns[stateid]+=1;
  file.clear();
  file.seekg(0);
  for(int j=0; j<ns[stateid]; j++)
    {
      for(int i=0; i<12; i++)
	file >> nums[i];
      datasigma[0].push_back(atof(nums[2])/mass[0]);
      datasigma[1].push_back(atof(nums[5])*mass[0]/(5.961*34.46e-4*1000));
      datasigma[2].push_back(sqrt(atof(nums[6])*atof(nums[6])+atof(nums[8])*atof(nums[8]))*mass[0]/(5.961*34.46e-4*1000));
      datasigma[3].push_back(stateid);
      datasigma[4].push_back(atof(nums[10])/atof(nums[11])-1);
      datasigma[5].push_back(atof(nums[3])/mass[0]);
      datasigma[6].push_back(atof(nums[4])/mass[0]);
      datasigma[7].push_back(atof(nums[0]));
      datasigma[8].push_back(atof(nums[1]));
      datasigma[9].push_back(7000.);
    }
  file.close();
  stateid++;
  
  //the ATLAS chic data must have the systematic uncertainty quadratically corrected for the luminosity uncertainty (included as a nuisance parameter)
  //chic2 ATLAS at 7 TeV
  file.open("data/ATLAS_chic2_csection.txt");
  if(file.is_open())
    while(getline(file,data))
      ns[stateid]+=1;
  file.clear();
  file.seekg(0);
  for(int j=0; j<ns[stateid]; j++)
    {
      for(int i=0; i<12; i++)
	file >> nums[i];
      aux=atof(nums[8])*atof(nums[8])-1.8e-2*atof(nums[5])*1.8e-2*atof(nums[5]);
      datasigma[0].push_back(atof(nums[2])/mass[1]);
      datasigma[1].push_back(atof(nums[5])*mass[1]/(1.5*5.961*19.2e-4*1000));
      datasigma[2].push_back(sqrt(atof(nums[6])*atof(nums[6])+aux)*mass[1]/(1.5*5.961*19.2e-4*1000));
      datasigma[3].push_back(stateid);
      datasigma[4].push_back(atof(nums[10])/atof(nums[11])-1);
      datasigma[5].push_back(atof(nums[3])/mass[1]);
      datasigma[6].push_back(atof(nums[4])/mass[1]);
      datasigma[7].push_back(atof(nums[0]));
      datasigma[8].push_back(atof(nums[1]));
      datasigma[9].push_back(7000.);
    }
  file.close();
  stateid++;

  //chic1 ATLAS at 7 TeV
  file.open("data/ATLAS_chic1_csection.txt");
  if(file.is_open())
    while(getline(file,data))
      ns[stateid]+=1;
  file.clear();
  file.seekg(0);
  for(int j=0; j<ns[stateid]; j++)
    {
      for(int i=0; i<12; i++)
	file >> nums[i];
      aux=atof(nums[8])*atof(nums[8])-1.8e-2*atof(nums[5])*1.8e-2*atof(nums[5]);
      datasigma[0].push_back(atof(nums[2])/mass[2]);
      datasigma[1].push_back(atof(nums[5])*mass[2]/(1.5*5.961*33.9e-4*1000));
      datasigma[2].push_back(sqrt(atof(nums[6])*atof(nums[6])+aux)*mass[2]/(1.5*5.961*33.9e-4*1000));
      datasigma[3].push_back(stateid);
      datasigma[4].push_back(atof(nums[10])/atof(nums[11])-1);
      datasigma[5].push_back(atof(nums[3])/mass[2]);
      datasigma[6].push_back(atof(nums[4])/mass[2]);
      datasigma[7].push_back(atof(nums[0]));
      datasigma[8].push_back(atof(nums[1]));
      datasigma[9].push_back(7000.);
    }
  file.close();
  stateid++;

  //jpsi CMS at 7 TeV
  file.open("data/CMS_jpsi_csection.txt");
  if(file.is_open())
    while(getline(file,data))
      ns[stateid]+=1;
  file.clear();
  file.seekg(0);
  for(int j=0; j<ns[stateid]; j++)
    {
      for(int i=0; i<12; i++)
	file >> nums[i];
      datasigma[0].push_back(atof(nums[2])/mass[3]);
      datasigma[1].push_back(atof(nums[5])*mass[3]/(5.961e-2*1000));
      datasigma[2].push_back(sqrt(atof(nums[6])*atof(nums[6])+atof(nums[8])*atof(nums[8]))*mass[3]/(5.961e-2*1000));
      datasigma[3].push_back(stateid);
      datasigma[4].push_back(atof(nums[10])/atof(nums[11])-1);
      datasigma[5].push_back(atof(nums[3])/mass[3]);
      datasigma[6].push_back(atof(nums[4])/mass[3]);
      datasigma[7].push_back(atof(nums[0]));
      datasigma[8].push_back(atof(nums[1]));
      datasigma[9].push_back(7000.);
    }
  file.close();
  stateid++;

  //chic2/chic1 ratio CMS at 7 TeV
  file.open("data/CMS_chic_ratio.txt");
  if(file.is_open())
    while(getline(file,data))
      ns[stateid]+=1;
  file.clear();
  file.seekg(0);
  for(int j=0; j<ns[stateid]; j++)
    {
      for(int i=0; i<12; i++)
	file >> nums[i];
      datasigma[0].push_back(atof(nums[2])/mass[3]);
      datasigma[1].push_back(atof(nums[5])*33.9/19.2);
      datasigma[2].push_back(sqrt(atof(nums[6])*atof(nums[6])+atof(nums[8])*atof(nums[8]))*33.9/19.2);
      datasigma[3].push_back(stateid);
      datasigma[4].push_back(atof(nums[10])-1);
      datasigma[5].push_back(atof(nums[3])/mass[3]);
      datasigma[6].push_back(atof(nums[4])/mass[3]);
      datasigma[7].push_back(atof(nums[0]));
      datasigma[8].push_back(atof(nums[1]));
      datasigma[9].push_back(7000.);
    }
  file.close();
  stateid++;

  //chib2/chib1 ratio CMS
  file.open("data/CMS_chib_ratio.txt");
  if(file.is_open())
    while(getline(file,data))
      ns[stateid]+=1;
  file.clear();
  file.seekg(0);
  for(int j=0; j<ns[stateid]; j++)
    {
      for(int i=0; i<12; i++)
	file >> nums[i];
      datasigma[0].push_back(atof(nums[2])/mass[6]);
      datasigma[1].push_back(atof(nums[5]));
      datasigma[2].push_back(sqrt(atof(nums[6])*atof(nums[6])+atof(nums[8])*atof(nums[8])));
      datasigma[3].push_back(stateid);
      datasigma[4].push_back(atof(nums[10])-1);
      datasigma[5].push_back(atof(nums[3])/mass[6]);
      datasigma[6].push_back(atof(nums[4])/mass[6]);
      datasigma[7].push_back(atof(nums[0]));
      datasigma[8].push_back(atof(nums[1]));
      datasigma[9].push_back(7000.);
    }
  file.close();
  stateid++;

  //ups(3S) CMS
  file.open("data/CMS_ups3_csection.txt");
  if(file.is_open())
    while(getline(file,data))
      ns[stateid]+=1;
  file.clear();
  file.seekg(0);
  for(int j=0; j<ns[stateid]; j++)
    {
      for(int i=0; i<12; i++)
	file >> nums[i];
      datasigma[0].push_back(atof(nums[2])/mass[4]);
      datasigma[1].push_back(atof(nums[5])*mass[4]/(2.4*2.18e4));
      datasigma[2].push_back(sqrt(atof(nums[6])*atof(nums[6])+atof(nums[8])*atof(nums[9]))*mass[4]/(2.4*2.18e4));
      datasigma[3].push_back(stateid);
      datasigma[4].push_back(abs(atof(nums[10])/atof(nums[11]))-1);
      datasigma[5].push_back(atof(nums[3])/mass[4]);
      datasigma[6].push_back(atof(nums[4])/mass[4]);
      datasigma[7].push_back(atof(nums[0]));
      datasigma[8].push_back(atof(nums[1]));
      datasigma[9].push_back(7000.);
    }
  file.close();
  stateid++;

  //ups(3S) ATLAS
  file.open("data/ATLAS_ups3_csection.txt");
  if(file.is_open())
    while(getline(file,data))
      ns[stateid]+=1;
  file.clear();
  file.seekg(0);
  for(int j=0; j<ns[stateid]; j++)
    {
      for(int i=0; i<12; i++)
	file >> nums[i];
      datasigma[0].push_back(atof(nums[2])/mass[4]);
      datasigma[1].push_back(atof(nums[5])*mass[4]/(2.18e4));
      datasigma[2].push_back(sqrt(atof(nums[6])*atof(nums[6])+atof(nums[8])*atof(nums[8]))*mass[4]/(2.18e4));
      datasigma[3].push_back(stateid);
      datasigma[4].push_back(abs(atof(nums[10])/atof(nums[11]))-1);
      datasigma[5].push_back(atof(nums[3])/mass[4]);
      datasigma[6].push_back(atof(nums[4])/mass[4]);
      datasigma[7].push_back(atof(nums[0]));
      datasigma[8].push_back(atof(nums[1]));
      datasigma[9].push_back(7000.);
    }
  file.close();
  stateid++;

  //ups(2S) CMS
  file.open("data/CMS_ups2_csection.txt");
  if(file.is_open())
    while(getline(file,data))
      ns[stateid]+=1;
  file.clear();
  file.seekg(0);
  for(int j=0; j<ns[stateid]; j++)
    {
      for(int i=0; i<12; i++)
	file >> nums[i];
      datasigma[0].push_back(atof(nums[2])/mass[5]);
      datasigma[1].push_back(atof(nums[5])*mass[5]/(2.4*1.93e4));
      datasigma[2].push_back(sqrt(atof(nums[6])*atof(nums[6])+atof(nums[8])*atof(nums[9]))*mass[5]/(2.4*1.93e4));
      datasigma[3].push_back(stateid);
      datasigma[4].push_back(abs(atof(nums[10])/atof(nums[11]))-1);
      datasigma[5].push_back(atof(nums[3])/mass[5]);
      datasigma[6].push_back(atof(nums[4])/mass[5]);
      datasigma[7].push_back(atof(nums[0]));
      datasigma[8].push_back(atof(nums[1]));
      datasigma[9].push_back(7000.);
    }
  file.close();
  stateid++;

  //ups(2S) ATLAS
  file.open("data/ATLAS_ups2_csection.txt");
  if(file.is_open())
    while(getline(file,data))
      ns[stateid]+=1;
  file.clear();
  file.seekg(0);
  for(int j=0; j<ns[stateid]; j++)
    {
      for(int i=0; i<12; i++)
	file >> nums[i];
      datasigma[0].push_back(atof(nums[2])/mass[5]);
      datasigma[1].push_back(atof(nums[5])*mass[5]/(1.93e4));
      datasigma[2].push_back(sqrt(atof(nums[6])*atof(nums[6])+atof(nums[8])*atof(nums[8]))*mass[5]/(1.93e4));
      datasigma[3].push_back(stateid);
      datasigma[4].push_back(abs(atof(nums[10])/atof(nums[11]))-1);
      datasigma[5].push_back(atof(nums[3])/mass[5]);
      datasigma[6].push_back(atof(nums[4])/mass[5]);
      datasigma[7].push_back(atof(nums[0]));
      datasigma[8].push_back(atof(nums[1]));
      datasigma[9].push_back(7000.);
    }
  file.close();
  stateid++;

  //ups(1S) CMS
  file.open("data/CMS_ups1_csection.txt");
  if(file.is_open())
    while(getline(file,data))
      ns[stateid]+=1;
  file.clear();
  file.seekg(0);
  for(int j=0; j<ns[stateid]; j++)
    {
      for(int i=0; i<12; i++)
	file >> nums[i];
      datasigma[0].push_back(atof(nums[2])/mass[6]);
      datasigma[1].push_back(atof(nums[5])*mass[6]/(2.4*2.48e4));
      datasigma[2].push_back(sqrt(atof(nums[6])*atof(nums[6])+atof(nums[8])*atof(nums[9]))*mass[6]/(2.4*2.48e4));
      datasigma[3].push_back(stateid);
      datasigma[4].push_back(abs(atof(nums[10])/atof(nums[11]))-1);
      datasigma[5].push_back(atof(nums[3])/mass[6]);
      datasigma[6].push_back(atof(nums[4])/mass[6]);
      datasigma[7].push_back(atof(nums[0]));
      datasigma[8].push_back(atof(nums[1]));
      datasigma[9].push_back(7000.);
    }
  file.close();
  stateid++;

  //ups(1S) ATLAS
  file.open("data/ATLAS_ups1_csection.txt");
  if(file.is_open())
    while(getline(file,data))
      ns[stateid]+=1;
  file.clear();
  file.seekg(0);
  for(int j=0; j<ns[stateid]; j++)
    {
      for(int i=0; i<12; i++)
	file >> nums[i];
      datasigma[0].push_back(atof(nums[2])/mass[6]);
      datasigma[1].push_back(atof(nums[5])*mass[6]/(2.48e4));
      datasigma[2].push_back(sqrt(atof(nums[6])*atof(nums[6])+atof(nums[8])*atof(nums[8]))*mass[6]/(2.48e4));
      datasigma[3].push_back(stateid);
      datasigma[4].push_back(abs(atof(nums[10])/atof(nums[11]))-1);
      datasigma[5].push_back(atof(nums[3])/mass[6]);
      datasigma[6].push_back(atof(nums[4])/mass[6]);
      datasigma[7].push_back(atof(nums[0]));
      datasigma[8].push_back(atof(nums[1]));
      datasigma[9].push_back(7000.);
    }
  file.close();
  stateid++;

  //13 TeV data
  //psiprime CMS
  file.open("data/CMS_psip_csection_13.txt");
  if(file.is_open())
    while(getline(file,data))
      ns[stateid]+=1;
  file.clear();
  file.seekg(0);
  for(int j=0; j<ns[stateid]; j++)
    {
      for(int i=0; i<12; i++)
	file >> nums[i];
      datasigma[0].push_back(atof(nums[2])/mass[0]);
      datasigma[1].push_back(atof(nums[5])*mass[0]/(7.9e-3*1000));
      datasigma[2].push_back(sqrt(atof(nums[6])*atof(nums[6])*atof(nums[5])*atof(nums[5])+atof(nums[8])*atof(nums[8])*atof(nums[5])*atof(nums[5]))*mass[0]/(100*7.9e-3*1000));
      datasigma[3].push_back(stateid);
      datasigma[4].push_back(atof(nums[10])/atof(nums[11])-1);
      datasigma[5].push_back(atof(nums[3])/mass[0]);
      datasigma[6].push_back(atof(nums[4])/mass[0]);
      datasigma[7].push_back(atof(nums[0]));
      datasigma[8].push_back(atof(nums[1]));
      datasigma[9].push_back(13000.);
    }
  file.close();
  stateid++;

  //jpsi CMS
  file.open("data/CMS_jpsi_csection_13.txt");
  if(file.is_open())
    while(getline(file,data))
      ns[stateid]+=1;
  file.clear();
  file.seekg(0);
  for(int j=0; j<ns[stateid]; j++)
    {
      for(int i=0; i<12; i++)
	file >> nums[i];
      datasigma[0].push_back(atof(nums[2])/mass[3]);
      datasigma[1].push_back(atof(nums[5])*mass[3]/(5.961e-2*1000));
      datasigma[2].push_back(sqrt(atof(nums[6])*atof(nums[6])*atof(nums[5])*atof(nums[5])+atof(nums[8])*atof(nums[8])*atof(nums[5])*atof(nums[5]))*mass[3]/(100*5.961e-2*1000));
      datasigma[3].push_back(stateid);
      datasigma[4].push_back(atof(nums[10])/atof(nums[11])-1);
      datasigma[5].push_back(atof(nums[3])/mass[3]);
      datasigma[6].push_back(atof(nums[4])/mass[3]);
      datasigma[7].push_back(atof(nums[0]));
      datasigma[8].push_back(atof(nums[1]));
      datasigma[9].push_back(13000.);
    }
  file.close();
  stateid++;

  //ups(3S) CMS
  file.open("data/CMS_ups3_csection_13.txt");
  if(file.is_open())
    while(getline(file,data))
      ns[stateid]+=1;
  file.clear();
  file.seekg(0);
  for(int j=0; j<ns[stateid]; j++)
    {
      for(int i=0; i<12; i++)
	file >> nums[i];
      datasigma[0].push_back(atof(nums[2])/mass[4]);
      datasigma[1].push_back(atof(nums[5])*mass[4]/(2.18e-2*1000));
      datasigma[2].push_back(sqrt(atof(nums[6])*atof(nums[6])*atof(nums[5])*atof(nums[5])+atof(nums[8])*atof(nums[9])*atof(nums[5])*atof(nums[5]))*mass[4]/(100*2.18e-2*1000));
      datasigma[3].push_back(stateid);
      datasigma[4].push_back(abs(atof(nums[10])/atof(nums[11]))-1);
      datasigma[5].push_back(atof(nums[3])/mass[4]);
      datasigma[6].push_back(atof(nums[4])/mass[4]);
      datasigma[7].push_back(atof(nums[0]));
      datasigma[8].push_back(atof(nums[1]));
      datasigma[9].push_back(13000.);
    }
  file.close();
  stateid++;

  //ups(2S) CMS
  file.open("data/CMS_ups2_csection_13.txt");
  if(file.is_open())
    while(getline(file,data))
      ns[stateid]+=1;
  file.clear();
  file.seekg(0);
  for(int j=0; j<ns[stateid]; j++)
    {
      for(int i=0; i<12; i++)
	file >> nums[i];
      datasigma[0].push_back(atof(nums[2])/mass[5]);
      datasigma[1].push_back(atof(nums[5])*mass[5]/(1.93e-2*1000));
      datasigma[2].push_back(sqrt(atof(nums[6])*atof(nums[6])*atof(nums[5])*atof(nums[5])+atof(nums[8])*atof(nums[9])*atof(nums[5])*atof(nums[5]))*mass[5]/(100*1.93e-2*1000));
      datasigma[3].push_back(stateid);
      datasigma[4].push_back(abs(atof(nums[10])/atof(nums[11]))-1);
      datasigma[5].push_back(atof(nums[3])/mass[5]);
      datasigma[6].push_back(atof(nums[4])/mass[5]);
      datasigma[7].push_back(atof(nums[0]));
      datasigma[8].push_back(atof(nums[1]));
      datasigma[9].push_back(13000.);
    }
  file.close();
  stateid++;

  //ups(1S) CMS
  file.open("data/CMS_ups1_csection_13.txt");
  if(file.is_open())
    while(getline(file,data))
      ns[stateid]+=1;
  file.clear();
  file.seekg(0);
  for(int j=0; j<ns[stateid]; j++)
    {
      for(int i=0; i<12; i++)
	file >> nums[i];
      datasigma[0].push_back(atof(nums[2])/mass[6]);
      datasigma[1].push_back(atof(nums[5])*mass[6]/(2.48e-2*1000));
      datasigma[2].push_back(sqrt(atof(nums[6])*atof(nums[6])*atof(nums[5])*atof(nums[5])+atof(nums[8])*atof(nums[9])*atof(nums[5])*atof(nums[5]))*mass[6]/(100*2.48e-2*1000));
      datasigma[3].push_back(stateid);
      datasigma[4].push_back(abs(atof(nums[10])/atof(nums[11]))-1);
      datasigma[5].push_back(atof(nums[3])/mass[6]);
      datasigma[6].push_back(atof(nums[4])/mass[6]);
      datasigma[7].push_back(atof(nums[0]));
      datasigma[8].push_back(atof(nums[1]));
      datasigma[9].push_back(13000.);
    }
  file.close();
}
