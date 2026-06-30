/// This is example of automatic TDC creation

#include "base/ProcMgr.h"
#include "hadaq/HldProcessor.h"
#include "hadaq/TdcProcessor.h"
#include "hadaq/TrbProcessor.h"


void after_create(hadaq::TrbProcessor* trb)
{
   printf("Called after all sub-components are created\n");

   if (!trb) return;
   trb->SetPrintErrors(10);

   printf("NUM TDC %d!\n", trb->NumSubProc());

   for (unsigned k=0;k<trb->NumSubProc();k++) {
      auto tdc = dynamic_cast<hadaq::TdcProcessor*>(trb->GetSubProc(k));
      if (!tdc) continue;

      printf("Configure %s!\n", tdc->GetName());

      // needed if tdc-ToT values are larger 60ns
      // arguments: delay, range (min, max)
      // DiRICH has a 20ns calibration pulse width
      // but it assumes a 312.5MHz Oscillator on it
      // The DiRICH5s1 has a 200MHz Calibration Oscillator
      // This gives a factor of 1.5625, so the calibration pulse has a length of
      // 31.25ns
      //tdc->SetToTRange(31.25, 40, 70);
      //tdc->SetToTRange(31.25, 30, 70);
      tdc->SetToTRange(20.345, 20, 90);
      //tdc->SetToTRange(30, 20, 90);
      // in new dirich5s design the correct PLL is included, so no correction necessary!

      // tdc->SetUseLastHit(true);
      for (unsigned nch=2;nch<tdc->NumChannels();nch++) {
        //tdc->SetRefChannel(nch, nch-1, 0xffffff, 2000,  -10., 10.); 
        //  tdc->SetRefChannel(nch, nch-1, 0x00001001, 2000,  -10., 10.);
      }
      //
      tdc->SetRefChannel(17, 18, 0x0079a2, 2000,  -10., 10.);
      tdc->SetRefChannel(18, 17, 0xffffff, 2000,  -10., 10.);

   }
}


void first()
{
  base::ProcMgr::instance()->SetCustomBinning("TDC_001001_Ch00_Tot", 600, 0, 300);
  
  // base::ProcMgr::instance()->SetRawAnalysis(true);
  base::ProcMgr::instance()->SetTriggeredAnalysis(true);

   // all new instances get this value
   // 1: component, 2: several per TDC, 3: several per channel 4: plus cross references
   base::ProcMgr::instance()->SetHistFilling(4);

   // this limits used for liner calibrations when nothing else is available
   hadaq::TdcMessage::SetFineLimits(31, 491);

   // default channel numbers and edges mask
   // 1 - use only rising edge, falling edge is ignore
   // 2   - falling edge enabled and fully independent from rising edge
   // 3   - falling edge enabled and uses calibration from rising edge
   // 4   - falling edge enabled and common statistic is used for calibration
      hadaq::TrbProcessor::SetDefaults(32, 2);

   // force usage of reftime from subevent header
   hadaq::TdcProcessor::SetTimeRefKind(3);

   // [min..max] range for TDC ids
   // hadaq::TrbProcessor::SetTDCRange(0x1000, 0x7FFF);

   // true indicates that dogma data are expected
   auto proc = new hadaq::TrbProcessor(0, nullptr, 5, true);

   //c27962
   //  fa798e
   //  b77982
   //  4e794e
   //proc->CreateTDC(0xc27962, 0xfa798e, 0xb77982, 0x4e794e, 0xca7990);

  // proc->CreateTDC(0x1001);                    // Dogs with MMCX, Trigger, Cherenkov 
   proc->CreateTDC(0x1048, 0x1049);
 //  proc->CreateTDC(0x1005, 0x1006);
    
 //  proc->CreateTDC(0x1003, 0x1004, 0x1005, 0x1006);    // Fiber1 X1&X2
 //  proc->CreateTDC(0x1007, 0x1008, 0x1009, 0x100a);    // Fiber1 Y1&Y2
 //  proc->CreateTDC(0x100b, 0x100c, 0x100d, 0x100e);    // Fiber2 X1&X2
 //  proc->CreateTDC(0x100f, 0x1010, 0x1011, 0x1012);   // Fiber2 Y1&Y2

 //  proc->CreateTDC(0x1013, 0x1014, 0x1015, 0x1016);   // DIRC Dogs
 //  proc->CreateTDC(0x1017, 0x1018);
 //  proc->CreateTDC(0x1019, 0x101a, 0x101b, 0x101d);
 //  proc->CreateTDC(0x101e, 0x101f, 0x1020, 0x1021);   
//   proc->CreateTDC(0x1022, 0x1023, 0x1024, 0x1025);
//   proc->CreateTDC(0x1026, 0x1027, 0x1028, 0x1029);
//   proc->CreateTDC(0x102a, 0x102b, 0x102c, 0x102d);
//   proc->CreateTDC(0x102e, 0x102f);                   // 28 Dogs for DIRC
   /*
   proc->CreateTDC(0x1030, 0x1031);   // 48 dogs
   proc->CreateTDC(0x1032, 0x1033, 0x1034, 0x1035);
   proc->CreateTDC(0x1036, 0x1037, 0x1038, 0x1039);   // 56 dogs
   */

//   proc->CreateTDC(0x103a, 0x103b);                   // Dogs for JStraw1
//   proc->CreateTDC(0x103c, 0x103d);                   // Dogs for JStraw2
//   proc->CreateTDC(0x103e, 0x103f);                   // Dogs for JStraw3
//   proc->CreateTDC(0x1040, 0x1041);                   // Dogs for KStraw1
//   proc->CreateTDC(0x1042, 0x1043);                   // Dogs for KStraw2
//   proc->CreateTDC(0x1044, 0x1045);                   // Dogs for KStraw3
   

   after_create(proc);

   // first parameter if filename  prefix for calibration files
   //     and calibration mode (empty string - no file I/O)
   // second parameter is hits count for autocalibration
   //     0 - only load calibration
   //    -1 - accumulate data and store calibrations only at the end
   //    -77 - accumulate data and store linear calibrations only at the end
   //    >0 - automatic calibration after N hits in each active channel
   //         if value ends with 77 like 10077 linear calibration will be calculated
   //    >1000000000 - automatic calibration after N hits only once, 1e9 excluding
   // third parameter is trigger type mask used for calibration
   //   (1 << 0xD) - special 0XD trigger with internal pulser, used also for TOT calibration
   //    0x3FFF - all kinds of trigger types will be used for calibration (excluding 0xE and 0xF)
   //   0x80000000 in mask enables usage of temperature correction
   proc->ConfigureCalibration("dogmacal_", 100000, 1 << 0xd);

   // only accept trigger type 0x1 when storing file
   // new hadaq::HldFilter(0x1);

   // create ROOT file store
   // base::ProcMgr::instance()->CreateStore("td.root");

   // 0 - disable store
   // 1 - std::vector<hadaq::TdcMessageExt> - includes original TDC message
   // 2 - std::vector<hadaq::MessageFloat>  - compact form, without channel 0, stamp as float (relative to ch0)
   // 3 - std::vector<hadaq::MessageDouble> - compact form, with channel 0, absolute time stamp as double
   base::ProcMgr::instance()->SetStoreKind(2);


   //base::ProcMgr::instance()->SetSecondName("second.C");

}

// extern "C" required by DABC to find function from compiled code


