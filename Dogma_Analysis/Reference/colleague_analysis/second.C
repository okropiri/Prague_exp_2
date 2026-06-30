#include <cstdio>
#include <fstream>
#include <iostream>
#include <vector>

#include "TTree.h"
#include "TH2D.h"
#include "TCanvas.h"

#include "TGo4Analysis.h"

#include "base/EventProc.h"
#include "base/Event.h"
#include "hadaq/TdcSubEvent.h"

class SecondProc : public base::EventProc {
protected:
   std::vector<std::string> fTdcIds;                     ///< TDC IDs for processing
   static constexpr unsigned int fNumChannels = 32;      ///< Number of channels per TDC
   static constexpr unsigned int fMaxTDCs = 2;           ///< Maximum number of TDCs
   static constexpr unsigned int fNumStrawDetectors = 3; ///< Number of straw detectors (X/Y pairs)
   static constexpr unsigned int fNumDircReadouts = 28;  ///< Number of DIRC readouts (16 or 64 fold PMTs)
   static constexpr unsigned int fScintTdcIndex = 6;     ///< Scintillator TDC index (0-based)
   static constexpr unsigned int fFiberTdcIndex = 7;     ///< Fiber TDC index (0-based)

   base::H1handle hNumHits;                        ///< Histogram with total number of hits
   base::H1handle hStart;                          ///< Scintillator start time
   
   base::H1handle hRF_ToT;                         ///< TOT RF Trig>
   base::H1handle hRF_T;                           ///< RF TDC >
                                                
   base::H1handle hScint[2];                       ///< Scint modules timing >
   base::H1handle hNcal[7];                        ///< Ncal modules timing >
   base::H1handle hVeto[2];                        ///< Veto timing >
   base::H1handle hGamma[2];                       ///< Gamma timing >

   base::H2handle hProfScint[2];                   // ToT vs T of small scint
   base::H2handle hProfNcal[8];                    // ToT vs T of Ncal
   base::H2handle hProfVeto[2];                    // ToT vs T of Veto
   base::H2handle hProfGamma[2];                   // ToT vs T of Gamma

   base::H2handle hProfLstilbene;                       // ToT vs T of large Stilbene
   base::H2handle hProfSstilbene;                       // ToT vs T of small Stilbene
    
   FILE* dataout;  ///< Output file for TDC data

public:
   SecondProc(const char* procname, const std::vector<std::string>& tdcids) : 
      base::EventProc(procname), fTdcIds(tdcids)
   {
      if (fTdcIds.size() != fMaxTDCs) {
         printf("Error: Expected %u TDCs, got %zu\n", fMaxTDCs, fTdcIds.size());
         fTdcIds.resize(fMaxTDCs, "");
      }
          // Replace the whole if-block with:
          if (fTdcIds.empty()) {
             printf("ERROR: no TDCs passed to SecondProc!\n");
             return; // or throw
          }
          // Do NOT resize!

      // Initialize histograms
      hNumHits = MakeH1("NumHits", "Number of hits", 100, 0, 500, "number");
      hRF_ToT = MakeH1("RF_ToT", "ToT of RF", 1000, 0, 100, "ns");
      hRF_T = MakeH1("RF_TDC", "TDC Time of RF", 1200,-1000, 200, "ns");

      for (unsigned c = 0; c < 2; c++){
         hScint[c] = MakeH1(TString::Format("Scint%d",c+1),
                            TString::Format("Scint_TDC_%d", c+1),
                           1000, -800, 200, "ns");
         hVeto[c] = MakeH1(TString::Format("Veto%d",c+1),
                           TString::Format("Veto_TDC_%d", c+1),
                           1000, -800, 200, "ns");
         hGamma[c] = MakeH1(TString::Format("Gamma%d",c+1),
                            TString::Format("Gamma_TDC%d",c+1),
                            1000, -800, 200, "ns");
      }

      for (unsigned c = 0; c < 7; c++){
         hNcal[c] = MakeH1(TString::Format("Ncal%d", c+1),
                           TString::Format("Ncal_TDC%d", c+1),
                           1000, -800, 200, "ns");

      }

      hProfScint[0] = MakeH2("ProfScint1", "Scint1_ToT_T", 1200, -1000, 200, 200, 0, 100, "T (ns); TOT (ns)");
      hProfScint[1] = MakeH2("ProfScint2", "Scint2_ToT_T", 1200, -1000, 200, 200, 0, 100, "T (ns); TOT (ns)");

      hProfVeto[0] = MakeH2("ProfVeto1", "Veto1_ToT_T", 1200, -1000, 200, 200, 0, 100, "T (ns); TOT (ns)");
      hProfVeto[1] = MakeH2("ProfVeto2", "Veto2_ToT_T", 1200, -1000, 200, 200, 0, 100, "T (ns); TOT (ns)");

      hProfGamma[0] = MakeH2("ProfGamma1", "Gamma1_ToT_T", 6000, -3000, 3000, 2000, 0, 1000, "T (ns); TOT (ns)");
      hProfGamma[1] = MakeH2("ProfGamma2", "Gamma2_ToT_T", 6000, -3000, 3000, 2000, 0, 1000, "T (ns); TOT (ns)");

      hProfNcal[0] = MakeH2("ProfNcal1", "Ncal1_ToT_T", 12000, -6000, 6000, 300, 0, 300, "T (ns); TOT (ns)");
      hProfNcal[1] = MakeH2("ProfNcal2", "Ncal2_ToT_T", 12000, -6000, 6000, 300, 0, 300, "T (ns); TOT (ns)");
      hProfNcal[2] = MakeH2("ProfNcal3", "Ncal3_ToT_T", 12000, -6000, 6000, 300, 0, 300, "T (ns); TOT (ns)");
      hProfNcal[3] = MakeH2("ProfNcal4", "Ncal4_ToT_T", 12000, -6000, 6000, 300, 0, 300, "T (ns); TOT (ns)");
      hProfNcal[4] = MakeH2("ProfNcal5", "Ncal5_ToT_T", 12000, -6000, 6000, 300, 0, 300, "T (ns); TOT (ns)");
      hProfNcal[5] = MakeH2("ProfNcal6", "Ncal6_ToT_T", 12000, -6000, 6000, 300, 0, 300, "T (ns); TOT (ns)");
      hProfNcal[6] = MakeH2("ProfNcal7", "Ncal7_ToT_T", 12000, -6000, 6000, 300, 0, 300, "T (ns); TOT (ns)");
      hProfNcal[7] = MakeH2("ProfNcal1_AttH", "Ncal1_ToT_T_AttH", 12000, -60000, 6000, 300, 0, 300, "T (ns); TOT (ns)");

      hProfLstilbene = MakeH2("ProfLstilbene", "Lstilbene_ToT_T", 12000, -6000, 6000, 300, 0, 300, "T (ns); TOT (ns)");
      hProfSstilbene = MakeH2("ProfSstilbene", "Sstilbene_ToT_T", 12000, -6000, 6000, 300, 0, 300, "T (ns); TOT (ns)");

      // Open output file 
      TGo4Analysis* analysis = TGo4Analysis::Instance();
      TString filename = analysis->GetInputFileName();
      TString filename1 = filename + ".dat";
      dataout = fopen(filename1.Data(), "a+");
      if (!dataout) {
         printf("Error: Cannot open output file %s\n", filename1.Data());
      }
   }

   virtual ~SecondProc()
   {
   //   if (dataout) fclose(dataout);
   }

   virtual void CreateBranch(TTree* t)
   {
      // Not used, as in original
   }

   virtual bool Process(base::Event* ev){
      static unsigned eventCounter = 0;
      eventCounter++;
       
      double num = 0;
      double times[fMaxTDCs][fNumChannels][2] = {{{0}}};          // times[tdc][ch][0=rising, 1=falling]
       
      // Process each TDC
      for (size_t t = 0; t < fTdcIds.size(); t++) {
         // <hadaq::TdcSubEventFloat*> returns relative time to Trigger
         // <hadaq::TdcSubEventDouble*> returns global time, wrap about 1 hour
         auto sub = dynamic_cast<hadaq::TdcSubEventFloat*>(ev->GetSubEvent(fTdcIds[t]));
         //   auto sub = dynamic_cast<hadaq::TdcSubEventDouble*>(ev->GetSubEvent(fTdcIds[t]));
           
         if (!sub) {
            // if (dataout) fprintf(dataout, "TDC %s not found\n", fTdcIds[t].c_str());
            continue;
         }
         // printf("t=%u, size %u, sub %p finding sub-events %p \n", t, fTdcIds.size(), sub, ev->GetSubEvent(fTdcIds[t]));
         // printf("TDC %s trigger time %8.9f\n Second",  fTdcIds[t].c_str(), sub->GetTriggerTime());
           
           
            if (dataout) {fprintf(dataout, "TDC %zu of total %zu TDCs sub %p finding sub-events %p TDC %s size %u\n",
            t, fTdcIds.size(),sub, ev->GetSubEvent(fTdcIds[t]), fTdcIds[t].c_str(), sub->Size());
            fprintf(dataout, "GlobalTriggerTime %8.9f\n", sub->GetTriggerTime());
            }
            
           
         for (unsigned cnt = 0; cnt < sub->Size(); cnt++) {
            auto& msg = sub->msg(cnt);
            float tm = msg.getStamp();
            unsigned chid = msg.getCh();
            bool isrising = msg.isRising();
            // cout<<"chid = "<<chid<<endl;
               
            num += 1;
               if (dataout) fprintf(dataout, "%zu %u %d %6.3f\n", t+1, chid, isrising, tm);
               
            if (chid < fNumChannels) {
               times[t][chid][isrising] = tm;
               // times[t][chid][isrising ? 0 : 1] = tm;
               // if(chid==1)printf("TDC_ID %u, Ch_ID %zu, Edge %d, Time %6.3f\n",t, chid, isrising, tm);
                   
            }
               
         }
      }
       
      if(eventCounter % 10 == 0){
         
        FillH1(hRF_T, times[1][0][1]);
        FillH1(hRF_ToT, times[1][0][0]-times[1][0][1]);

        // Fill the TDC time of NCal
        FillH1(hNcal[0], times[1][2][1]);
        FillH1(hNcal[1], times[1][3][1]);
        FillH1(hNcal[2], times[1][4][1]);
        FillH1(hNcal[3], times[1][5][1]);
        FillH1(hNcal[4], times[1][6][1]);
        FillH1(hNcal[5], times[1][7][1]);
        FillH1(hNcal[6], times[1][8][1]);

        //  Fill the TDC time of Veto detectors

        FillH1(hVeto[0], times[1][10][1]);
        FillH1(hVeto[1], times[1][11][1]);

        // Fill the TDC time of small scintillator

        FillH1(hScint[0], times[1][13][1]);
        FillH1(hScint[1], times[1][14][1]);

        // Fill the TDC time of Gamma detector

        FillH1(hGamma[0], times[1][16][1]);
        FillH1(hGamma[1], times[1][17][1]);
      

        // fill ToT vs T of Ncal
         if(times[1][2][0] != 0 && times[1][2][1] != 0){
            FillH2(hProfNcal[0], times[1][2][1], times[1][2][0]-times[1][2][1]);
         }
         if(times[1][3][0] != 0 && times[1][3][1] != 0){
            FillH2(hProfNcal[1], times[1][3][1], times[0][3][0]-times[1][3][1]);
         }
         if(times[1][4][0] != 0 && times[1][4][1] != 0){
            FillH2(hProfNcal[2], times[1][4][1], times[0][4][0]-times[1][4][1]);
         }
         if(times[1][5][0] != 0 && times[1][5][1] != 0){
            FillH2(hProfNcal[3], times[1][5][1], times[0][5][0]-times[1][5][1]);
         }
         if(times[1][6][0] != 0 && times[1][6][1] != 0){
            FillH2(hProfNcal[4], times[1][6][1], times[1][6][0]-times[1][6][1]);
         }
         
         if(times[1][7][0] != 0 && times[1][7][1] != 0){
            FillH2(hProfNcal[5], times[1][7][1], times[1][7][0]-times[1][7][1]);
         }
         if(times[1][8][0] != 0 && times[1][8][1] != 0){
            FillH2(hProfNcal[6], times[1][8][1], times[1][8][0]-times[1][8][1]);
         }

         if(times[1][9][0] != 0 && times[1][9][1] != 0){
            FillH2(hProfNcal[7], times[1][9][1], times[1][9][0]-times[1][9][1]);
         }
         
         // Fill ToT vs T of Veto
         if(times[1][11][1] != 0 && times[1][11][0] != 0){
           FillH2(hProfVeto[0], times[1][11][1], times[1][11][0]-times[1][11][1]);
         }

         if(times[1][12][1] != 0 && times[1][12][0] != 0){
           FillH2(hProfVeto[1], times[1][12][1], times[1][12][0]-times[1][12][1]);
         }

         // Fill ToT vs T of Scintillator
         if(times[1][14][1] != 0 && times[1][14][0] != 0){
           FillH2(hProfScint[0], times[1][14][1], times[1][14][0]-times[1][14][1]);
         }

         if(times[1][15][1] != 0 && times[1][15][0] != 0){
           FillH2(hProfVeto[1], times[1][15][1], times[1][15][0]-times[1][15][1]);
         }

         // Fill ToT vs T of Gamma
         if(times[1][17][1] != 0 && times[1][17][0] != 0){
           FillH2(hProfGamma[0], times[1][17][1], times[1][17][0]-times[1][17][1]);
         }

         if(times[1][18][1] != 0 && times[1][18][0] != 0){
           FillH2(hProfGamma[1], times[1][18][1], times[1][18][0]-times[1][18][1]);
         }

         if(times[1][20][1] != 0 && times[1][20][0] != 0){
           FillH2(hProfSstilbene, times[1][20][1], times[1][20][0]-times[1][20][1]);
         }
         if(times[1][22][1] != 0 && times[1][22][0] != 0){
           FillH2(hProfLstilbene, times[1][22][1], times[1][22][0]-times[1][22][1]);
         }
           
         FillH1(hNumHits, num);
         return true;
      }
   }
};

void second(){
   std::vector<std::string> tdcids = {
       "TDC_001048", "TDC_001049"
      };
      
   
   new SecondProc("custom", tdcids);
}
