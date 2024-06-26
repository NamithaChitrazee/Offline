//
// Make plots from the ntuple made by Analyses/src/stoppingTarget00_plugin.cc
//
//
// Original author Rob Kutschke
//

{

  // With this you can reinvoke the script without exiting root.
  gROOT->Reset();

  // Enable printing and prompt at end of each page
  bool doprint = true;
  bool prompt = false;

  // Which pages to do; the scatterplots are very big.
  bool page01 = true;
  bool page02 = true;   // Scatterplot
  bool page02a = true;
  bool page03 = true;
  bool page04 = true;
  bool page05 = true;
  bool page06 = true;   // Scatterplot
  bool page07 = true;   // Scatterplot
  bool page08 = true;
  bool page09 = true;
  bool page10 = true;

  // Number of entries in the true scatterplots.
  int nDrawPage2=50000;
  int nDrawPage2d=75000;
  int nDrawPage6=50000;  // For pages 6.
  int nDrawPage7=75000;  // For pages 7.

  // Get rid of grey background (ugly for print out).
  gROOT->SetStyle("Plain");

  // Statistics box for histograms should include all of:
  // number of Entries, Mean, Rms, Underflows, Overflows
  gStyle->SetOptStat("emruo");

  //gStyle->SetMarkerStyle(26);
  //gStyle->SetMarkerSize(0.625);
  gStyle->SetMarkerColor(kBlue);
  gStyle->SetHistLineColor(kBlue);

  gROOT->ForceStyle();
  //TString basename("stoppingTarget00_20000");
  //TString basename("stoppingTarget00");
  TString basename("stoppingTarget00_Full");

  // Open the input file that contains histograms and ntuples
  // made by ReadBack.cc
  TFile* file = new TFile( basename + ".root");

  TString psfile( basename + ".ps");

  TH1F* hStopFoil;     file->GetObject("stopping/hStopFoil",  hStopFoil);
  TH1F* hnSimPart;     file->GetObject("stopping/hnSimPart", hnSimPart);
  TH1F* hnSimPartZoom; file->GetObject("stopping/hnSimPartZoom",     hnSimPartZoom);

  TNtuple* nt; file->GetObject("stopping/nt",nt);

  // Some standard cuts
  TCut decayed  = "scode==14";
  TCut stopped  = "scode==32";
  TCut inTarget = "scode==32&&sfoil>-1";           // Stopped in the target
  TCut reachedEnd  = "scode==32&&sz>15800";
  TCut outOfTarget = "sfoil==-1";                  // Stopped out of the target
  TCut interactOutofTarget = stopped&&outOfTarget&&!reachedEnd; // Interacted out of the target

  // Some computed ntuple variables
  char * cr    = "sqrt(cx*cx+cy*cy)";     // Radius at entry point
  char * sr    = "sqrt(sx*sx+sy*sy)";     // Radius at stopping point
  char * rvd9  = "sqrt(x9*x9+y9*y9)";     // Radius at vd 9.
  char * rvd10 = "sqrt(x10*x10+y10*y10)"; // Radius at vd 10.

  // Open a new canvas on the screen.
  TCanvas *c = new TCanvas("c", "Plots from " + basename, 900, 900 );

  // Open a multi-page output postscript file .
  if ( doprint ) c->Print( psfile+"[");

  if ( page01 ){
    // Split the canvas into 4 pads.
    c->Divide(1,2);

    c->cd(1); hStopFoil->Draw("H9");
    c->cd(2); hnSimPart->Draw("H9");

    // Flush page to screen
    c->Modified();
    c->Update();

    // Add this canvas to the postscript file.
    if ( doprint )  c->Print(psfile);

    // Prompt and wait for response before continuing.
    if ( prompt) {
      cerr << "Double click in the last active pad to continue:" ;
      gPad->WaitPrimitive();
      cerr << endl;
    }

    // Clear canvas in preparation for next page
    gPad->Clear();
    c->cd(0);
    c->Clear();
  }

  if ( page02 ){

    c->Divide(2,2);
    c->cd(1) ;
    TH1F* frame = c_1->DrawFrame(-250.,-250.,250.,250.);
    nt->Draw( "cx:cy","","PSAME",nDrawPage2);
    frame->SetTitle("y vs x at Entry;(mm);(mm)");

    c->cd(2);
    frame = c_2->DrawFrame(-250.,-250.,250.,250.);
    nt->Draw( "x9:y9","","PSAME",nDrawPage2);
    frame->SetTitle("y vs x at VD 9;(mm);(mm)");

    c->cd(3);
    frame = c_3->DrawFrame(-250.,-250.,250.,250.);
    nt->Draw( "x10:y10","","PSAME",nDrawPage2);
    frame->SetTitle("y vs x at VD 10;(mm);(mm)");

    c->cd(4);
    frame = c_4->DrawFrame(-250.,-250.,250.,250.);
    nt->Draw( "sx:sy",inTarget,"PSAME",nDrawPage2d);
    frame->SetTitle("y vs x at Stopping Point in Target;(mm);(mm)");

    // Flush page to screen
    c->Modified();
    c->Update();

    if ( doprint ) c->Print(psfile);

    // Prompt and wait for response before continuing.
    if ( prompt ){
      cerr << "Double click in the last active pad to continue:" ;
      gPad->WaitPrimitive();
      cerr << endl;
    }

    // Clear canvas in preparation for next page.
    gPad->Clear();
    c->cd(0);
    c->Clear();
    c->Modified();
    c->Update();
  }

  if ( page02a ){

    c->Divide(2,2);
    c->cd(1) ;
    TH2F* hxyCreation = new TH2F( "hxyCreation",
                                  "y vs x at Creation All Tracks;(mm);(mm)",
                                  50, -250., 250.,
                                  50, -250., 250. );
    nt->Project( "hxyCreation", "cx:cy" );
    hxyCreation->Draw("Box");

    c->cd(2) ;
    TH2F* hxyvd9 = new TH2F( "hxyvd9",
                                  "y vs x at Virtual Detector 9 All Tracks;(mm);(mm)",
                                  50, -250., 250.,
                                  50, -250., 250. );
    nt->Project( "hxyvd9", "x9:y9", "sz>5499.89" );
    hxyvd9->Draw("Box");

    c->cd(3) ;
    TH2F* hxyvd10 = new TH2F( "hxyvd10",
                                  "y vs x at Virtual Detector 10 All Tracks;(mm);(mm)",
                                  50, -250., 250.,
                                  50, -250., 250. );
    nt->Project( "hxyvd10", "x10:y10", "sz>6301.11" );
    hxyvd10->Draw("Box");

    c->cd(4) ;
    TH2F* hxyTarget = new TH2F( "hxyTarget",
                                  "y vs x at Stopping Point in Target;(mm);(mm)",
                                  50, -250., 250.,
                                  50, -250., 250. );
    nt->Project( "hxyTarget", "sx:sy", inTarget );
    hxyTarget->Draw("Box");

    // Flush page to screen
    c->Modified();
    c->Update();

    if ( doprint ) c->Print(psfile);

    // Prompt and wait for response before continuing.
    if ( prompt ){
      cerr << "Double click in the last active pad to continue:" ;
      gPad->WaitPrimitive();
      cerr << endl;
    }

    // Clear canvas in preparation for next page.
    gPad->Clear();
    c->cd(0);
    c->Clear();
    c->Modified();
    c->Update();
  }

  if ( page03 ) {

    TH1F* hr0          = new TH1F( "hr0",
                                   "Radius at entry", 100, 0., 200. );
    TH1F* hr0Stopped   = new TH1F( "hr0Stopped",
                                   "Radius at entry iff stopped in Foil", 100, 0., 200. );
    TH1F* hrvd9        = new TH1F( "hrvd9",
                                   "Radius at VD 9", 100, 0., 400. );
    TH1F* hrvd9Stopped = new TH1F( "hrvd9Stopped",
                                   "Radius at VD 9 iff stopped in Foil", 100, 0., 400. );
    TH1F* hrvd10       = new TH1F( "hrvd10",
                                   "Radius at VD 10", 100, 0., 400. );

    nt->Project( "hr0",          cr );
    nt->Project( "hr0Stopped",   cr, inTarget );
    nt->Project( "hrvd9",        rvd9, "sz>5200" );
    nt->Project( "hrvd9Stopped", rvd9, inTarget );
    nt->Project( "hrvd10",       rvd10, "sz>6300" );

    double max1 = hrvd9->GetMaximum()*1.05;
    hr0->SetMaximum(max1);
    hr0Stopped->SetMaximum(max1);
    hrvd9->SetMaximum(max1);
    hrvd9Stopped->SetMaximum(max1);
    hrvd10->SetMaximum(max1);

    hr0Stopped->SetLineColor(kRed);
    hrvd9Stopped->SetLineColor(kRed);

    c->Divide(2,3);

    c->cd(1); hr0->Draw("H9");
              hr0Stopped->Draw("H9SAME");
    c->cd(2); hr0Stopped->Draw("H9SAME");
    c->cd(3); hrvd9->Draw("H9");
               hrvd9Stopped->Draw("H9SAME");
    c->cd(4); hrvd9Stopped->Draw("H9");
    c->cd(5); hrvd10->Draw("H9");

    // Flush page to screen
    c->Modified();
    c->Update();
    if ( prompt ){
      cerr << "Double click in the last active pad to continue:" ;
      gPad->WaitPrimitive();
      cerr << endl;
    }

    if ( doprint ) c->Print(psfile);

    // Clear canvas in preparation for next page
    gPad->Clear();
    c->cd(0);
    c->Clear();
  }

  if ( page04 ) {

    TH1F* hp0          = new TH1F( "hp0",
                                   "Momentum at entry", 100, 0., 120. );
    TH1F* hp0Stopped   = new TH1F( "hp0Stopped",
                                   "Momentum at entry iff stopped in Foil", 100, 0., 120. );

    TH1F* hpt0          = new TH1F( "hpt0",
                                   "pT at entry", 100, 0., 60. );
    TH1F* hpt0Stopped   = new TH1F( "hpt0Stopped",
                                   "pT at entry iff stopped in Foil", 100, 0., 60. );

    nt->Project( "hp0",         "cp" );
    nt->Project( "hp0Stopped",  "cp", inTarget );
    nt->Project( "hpt0",        "cpt" );
    nt->Project( "hpt0Stopped", "cpt", inTarget );

    hp0Stopped->SetLineColor(kRed);
    hpt0Stopped->SetLineColor(kRed);

    c->Divide(2,2);

    c->cd(1); hp0->Draw("H9");
              hp0Stopped->Draw("H9SAME");
    c->cd(2); hp0Stopped->Draw("H9SAME");

    c->cd(3); hpt0->Draw("H9");
              hpt0Stopped->Draw("H9SAME");
    c->cd(4); hpt0Stopped->Draw("H9SAME");

    // Flush page to screen
    c->Modified();
    c->Update();
    if ( prompt ){
      cerr << "Double click in the last active pad to continue:" ;
      gPad->WaitPrimitive();
      cerr << endl;
    }

    // Clear canvas in preparation for next page
    gPad->Clear();
    c->cd(0);
    c->Clear();
  }



  if ( page05 ) {
    TH1F* hszDecay = new TH1F( "hszDecay", "z at Decay", 300, 0., 16000. );
    TH1F* hszFoils = new TH1F( "hszFoils", "z at Stopping Point in Foil", 120, 5400., 6400. );
    TH1F* hszOther = new TH1F( "hszOther", "z at Interaction out of Foil", 300, -4000., 16000. );
    TH1F* hszEnd   = new TH1F( "hszEnd",   "Reached end of DS", 300, 15000., 16000. );
    nt->Project( "hszDecay",  "sz", decayed );
    nt->Project( "hszFoils",  "sz", inTarget);
    nt->Project( "hszOther",  "sz", interactOutofTarget );
    nt->Project( "hszEnd",    "sz", reachedEnd );

    c->Divide(2,2);

    hszDecay->GetXaxis()->SetNdivisions(504);
    hszFoils->GetXaxis()->SetNdivisions(504);
    hszOther->GetXaxis()->SetNdivisions(504);
    hszEnd->GetXaxis()->SetNdivisions(504);

    c->cd(1); hszDecay->Draw("H9");
    c->cd(2); hszFoils->Draw("H9");
    c->cd(3); hszOther->Draw("H9");
    c->cd(4); hszEnd->Draw("H9");

    // Flush page to screen
    c->Modified();
    c->Update();

    if ( doprint ) c->Print(psfile);

    // Prompt and wait for response before continuing.
    if ( prompt ){
      cerr << "Double click in the last active pad to continue:" ;
      gPad->WaitPrimitive();
      cerr << endl;
    }

    // Clear canvas in preparation for next page.
    gPad->Clear();
    c->cd(0);
    c->Clear();
    c->Modified();
    c->Update();

  }

  if ( page06 ){

    c->Divide(3,3);

    for ( int ifoil=0; ifoil<9; ++ifoil ){
      int ipad = ifoil +1;

      ostringstream cut;
      cut << "sfoil==" << ifoil <<"&&scode==32";

      ostringstream title;
      title << "y vs x If Stopped in Foil " << ifoil << ";(mm);(mm)";

      c->cd(ipad);
      TH1F* frame = gPad->DrawFrame(-110.,-110.,110.,110.);
      nt->Draw( "sx:sy", cut.str().c_str(), "PSAME", nDrawPage6 );
      frame->SetTitle( title.str().c_str() );

    }

    // Flush page to screen
    c->Modified();
    c->Update();

    if ( doprint ) c->Print(psfile);

    // Prompt and wait for response before continuing.
    if ( prompt ){
      cerr << "Double click in the last active pad to continue:" ;
      gPad->WaitPrimitive();
      cerr << endl;
    }

    // Clear canvas in preparation for next page.
    gPad->Clear();
    c->cd(0);
    c->Clear();
    c->Modified();
    c->Update();
  }

  if ( page07 ){

    c->Divide(3,3);

    for ( int ifoil=9; ifoil<17; ++ifoil ){
      int ipad = ifoil-8;

      ostringstream cut;
      cut << "sfoil==" << ifoil <<"&&scode==32";

      ostringstream title;
      title << "y vs x If Stopped in Foil " << ifoil << ";(mm);(mm)";

      c->cd(ipad);
      TH1F* frame = gPad->DrawFrame(-110.,-110.,110.,110.);
      nt->Draw( "sx:sy", cut.str().c_str(), "PSAME", nDrawPage7 );
      frame->SetTitle( title.str().c_str() );

    }

    // Flush page to screen
    c->Modified();
    c->Update();

    if ( doprint ) c->Print(psfile);

    // Prompt and wait for response before continuing.
    if ( prompt ){
      cerr << "Double click in the last active pad to continue:" ;
      gPad->WaitPrimitive();
      cerr << endl;
    }

    // Clear canvas in preparation for next page.
    gPad->Clear();
    c->cd(0);
    c->Clear();
    c->Modified();
    c->Update();
  }


  if ( page08 ){

    c->Divide(3,3);

    for ( int ifoil=0; ifoil<9; ++ifoil ){
      int ipad = ifoil +1;

      ostringstream name;
      name << "foilxy_" << ifoil;

      ostringstream cut;
      cut << "sfoil==" << ifoil <<"&&scode==32";

      ostringstream title;
      title << "y vs x If Stopped in Foil " << ifoil << ";(mm);(mm)";

      c->cd(ipad);
      TH2F* hist = new TH2F( name.str().c_str(),
                             title.str().c_str(),
                             30, -120., 120.,
                             30, -120., 120. );

      nt->Project( name.str().c_str(), "sx:sy", cut.str().c_str() );

      hist->Draw("BOX");

    }

    // Flush page to screen
    c->Modified();
    c->Update();

    if ( doprint ) c->Print(psfile);

    // Prompt and wait for response before continuing.
    if ( prompt ){
      cerr << "Double click in the last active pad to continue:" ;
      gPad->WaitPrimitive();
      cerr << endl;
    }

    // Clear canvas in preparation for next page.
    gPad->Clear();
    c->cd(0);
    c->Clear();
    c->Modified();
    c->Update();
  }

  if ( page09 ){

    c->Divide(3,3);

    for ( int ifoil=9; ifoil<17; ++ifoil ){
      int ipad = ifoil-8;

      ostringstream name;
      name << "foilxy_" << ifoil;

      ostringstream cut;
      cut << "sfoil==" << ifoil <<"&&scode==32";

      ostringstream title;
      title << "y vs x If Stopped in Foil " << ifoil << ";(mm);(mm)";

      c->cd(ipad);
      TH2F* hist = new TH2F( name.str().c_str(),
                             title.str().c_str(),
                             30, -120., 120.,
                             30, -120., 120. );

      nt->Project( name.str().c_str(), "sx:sy", cut.str().c_str() );

      hist->Draw("BOX");

    }

    // Flush page to screen
    c->Modified();
    c->Update();

    if ( doprint ) c->Print(psfile);

    // Prompt and wait for response before continuing.
    if ( prompt ){
      cerr << "Double click in the last active pad to continue:" ;
      gPad->WaitPrimitive();
      cerr << endl;
    }

    // Clear canvas in preparation for next page.
    gPad->Clear();
    c->cd(0);
    c->Clear();
    c->Modified();
    c->Update();
  }

  if ( page10 ){

    c->Divide(2,2);

    c->cd(1);
    TH2F* hFoils1 = new TH2F( "hnFoilsVssFoil",
                           "Number of hit foils vs stopping foil;Stopping Foil Number; Number of Hit Foils",
                           22, 0., 22.,
                           22, 0., 22. );

    nt->Project( "hnFoilsVssFoil", "nfoils:sfoil", "sfoil>-1&&scode==32" );
    hFoils1->Draw("BOX");

    c->cd(2);
    TH1F* hFoils2 = new TH1F( "hFoils2",
                              "Number of hit foils for particles stopping elsewhere; Foil Number",
                              22, -1., 21. );
    nt->Project( "hFoils2", "nfoils","scode==32&&sfoil==-1");
    hFoils2->Draw("H9");

    c->cd(3);
    TH1F* hFoils3 = new TH1F( "hFoils3",
                              "Number of hit foils for particles that decay; Foil Number",
                              22, -1., 21. );
    nt->Project( "hFoils3", "nfoils","scode==14&&sfoil==-1");
    hFoils3->Draw("H9");

    // Flush page to screen
    c->Modified();
    c->Update();

    if ( doprint ) c->Print(psfile);

    // Prompt and wait for response before continuing.
    if ( prompt ){
      cerr << "Double click in the last active pad to continue:" ;
      gPad->WaitPrimitive();
      cerr << endl;
    }

    // Clear canvas in preparation for next page.
    gPad->Clear();
    c->cd(0);
    c->Clear();
    c->Modified();
    c->Update();
  }

  // Close the postscript file.
  if ( doprint ) c->Print(psfile+"]");

}
