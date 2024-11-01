double PulseShape(double *xx,double *par)
{
        double &x=*xx, &a0=par[0], &am=par[1], &tm=par[2], &tau=par[3], &p=par[4];

        if( x<=tm-tau*log(p+1) ) return a0;

        double val, e;
        e=exp((tm-x)/tau);
        val=a0 + am*e*pow((1+p-e)/p,p);;
        if( val>1023 ) return 1023;
        return val;
}

void get_npe(TString list)
{	
        TString inDir = "/home/ovtin/development/TB_SND_ASHIPH/TBsummer2024/";
        TString outDir = "/home/ovtin/development/TB_SND_ASHIPH/TBsummer2024/results/";
	TString outfile = list + "_res.root";
	bool TimeDist = 0;
	bool waveformdraw = 1;

	int Nevexbeamdata, Nevtrack, Nevdaq, NevT;
	TString NEXBEAMDATA, NDAQ;
	TString inDirData; 
	TString intrackData;
	TString inDAQData;
	TString inTData;

	TChain *exbeamdata = new TChain("exbeamdata");
	TChain *track = new TChain("track");
	TChain *daq = new TChain("daq");
	TChain *T = new TChain("T");

	string line;	
	ifstream fin(inDir + "list_of_runs/" + list);
        if(fin.is_open())
	{
		while (getline(fin, line))
    		{
			if( line.size()>=1 && line[0]!='#' )
			{
				std::istringstream i_str(line);
				i_str>>NEXBEAMDATA>>NDAQ;
				
				inDirData = inDir + NEXBEAMDATA + ".root" ; 
				intrackData = inDir + "GEM/track_" + NDAQ + ".root";
				inDAQData = inDir + "GEM/" + NDAQ + ".root";                    
				if(TimeDist== 1) inTData = inDir + "TIME/" + NEXBEAMDATA +"_time.root";   

				exbeamdata->Add(inDirData);
				track->Add(intrackData);
				daq->Add(inDAQData);
				if(TimeDist== 1) T->Add(inTData);

				Nevexbeamdata=daq->GetEntries();
				Nevtrack=track->GetEntries();
				Nevdaq=daq->GetEntries();
				if(TimeDist== 1) NevT=T->GetEntries();

				if(TimeDist== 1 && Nevexbeamdata==Nevtrack && Nevexbeamdata==Nevdaq && Nevexbeamdata==NevT)
				{
					
					cout<<"Number of events is good"<<endl;
				}
				else if(TimeDist==0 && Nevexbeamdata==Nevtrack && Nevexbeamdata==Nevdaq)
				{
					cout<<"Number of events is good"<<endl;	
				}
				else
				{
					cout<<"Number of events is bad"<<endl;
					break;
				}
   			}
		}
	}
	else{
		cout << "Can't open input text file !"<<endl;
	}
   	fin.close();	
	
	exbeamdata->Print();    
	track->Print();
	daq->Print();
	if(TimeDist== 1) T->Print();

	exbeamdata->AddFriend("track");
	exbeamdata->AddFriend("daq");
	if(TimeDist== 1) exbeamdata->AddFriend("T");
	
	//----------------------------------------------------------------------
	auto f = TFile::Open(outDir + outfile,"RECREATE");        
	
	TProfile* ampped = new TProfile("ampped", "a(t) for Cnt (pedestals) ", 1024, 0. , 1024.);
	TProfile* amp0 = new TProfile("amp0", "a(t) for Cnt", 1024, 0., 1024.);
		
	TH2F* htimeAmp = new TH2F("htimeAmp", "ha(t) for Cnt",1024,0,1024,4000,0,4000);
	TH2F* htimeAmpSelected = new TH2F("htimeAmpSelected", "a(t) for Cnt",1024,0,1024,5000,-1000,4000);

	TH2F* hxyGEMs0 = new TH2F("hxyGEMs0", "GEMs map of XY", 120,0,120,80,0,80);
	TH2F* hxyGEMs = new TH2F("hxyGEMs", "GEMs map of XY", 120,0,120,80,0,80);
	TH2F* hxyGEMswoXYcuts = new TH2F("hxyGEMswoXYcuts", "GEM map of XY wo cuts on XY", 120,0,120,80,0,80);
	TH2F* hxyGEMsXYWLS = new TH2F("hxyGEMsXYWLS", "GEM map of XY wo cuts on XY", 120,0,120,80,0,80);

	TH1F* hchtrg1 = new TH1F("hchtrg1", "Trg1 amplitude", 4096,0,4096);
	TH1F* hchtrg2 = new TH1F("hchtrg2", "Trg2 amplitude", 4096,0,4096);

	TH1F* hAmpSig = new TH1F("hAmpSig", "Counter amplitude", 2000,-100,1900);
	TH1F* hAmpPed = new TH1F("hAmpPed", "Counter amplitude", 2000,-100,1900);
        TH1F* hAmpSigInt = new TH1F("hAmpSigInt", "Counter integral", 1000,0,50000);
        TH1F* hAmpPedInt = new TH1F("hAmpPedInt", "Counter integral", 1000,0,50000);
	
        TH1F* hnpeAmpSig = new TH1F("hnpeAmpSig", "Nph.e. in counter", 400,0,200);
	TH1F* hnpeIntSig = new TH1F("hnpeIntSig", "Nph.e. in counter", 400,0,200);

	TH1F* hcal = new TH1F("hcal", "Energy of calorimeter", 40000,-150,1500);
	TH1F* hcalCuts = new TH1F("hcalCuts", "Energy of calorimeter with cuts", 40000,-150,1500);

	TH1F* hdtTrg = new TH1F("hdtTrg", "Time distribution between triggers 1 and 2", 500,-50.0,50.0);
	TH1F* hdtCntTrg1 = new TH1F("hdtCntTrg1", "Time distribution between Cnt and trigger 1", 500,-50.0,50.0);
	TH1F* hdtCntTrg2 = new TH1F("hdtCntTrg2", "Time distribution between Cnt and trigger 2", 500,-50.0,50.0);
	auto hdtCntTrg1AmpCnt = new TProfile("hdtCntTrg1AmpCnt","#Delta T Cnt and trigger 1 VS amplitude Cnt",100,0,500);
	auto hdtCntTrg2AmpCnt = new TProfile("hdtCntTrg2AmpCnt","#Delta T Cnt and trigger 2 VS amplitudes Cnt",100,0,500);

	//------------------------------------------------------------------------

	TString DATA;
	TString CONDITION;

	Double_t X1 = 61.0, X2 = 68.0; 
	Double_t Y1 = 27.0, Y2 = 37.0; //mm Geometrical borders
	
	Int_t t1L=200, t1R=10000, t2L=200, t2R=10000; //trig amps borders

	Int_t CEnL=500, CEnR=10000; //calorimeter energy borders

	TString by1;
	by1.Form("%g", Y1);
	TString by2;
	by2.Form("%g", Y2);
	TString bx1;
	bx1.Form("%g", X1);
	TString bx2;
	bx2.Form("%g", X2);
	TString bt1L;
	bt1L.Form("%d", t1L);
	TString bt1R;
	bt1R.Form("%d", t1R);
	TString bt2L;
	bt2L.Form("%d", t2L);
	TString bt2R;
	bt2R.Form("%d", t2R);
	TString bCEnL;
	bCEnL.Form("%d", CEnL);
	TString bCEnR;
	bCEnR.Form("%d", CEnR);
	
	//-----------------------------------------------------------------
	DATA = "lecroy2249.ch0>>hcal";
	TString CONDITIONwoCutsOnEn = "track.ys>"+by1+"   &&track.ys<"+by2+"  &&track.xs>"+bx1+"   &&track.xs<"+bx2+"  &&(chtrg1r.ped-chtrg1r.min)>"+bt1L+"&&(chtrg1r.ped-chtrg1r.min)<"+bt1R+"&&(chtrg2r.ped-chtrg2r.min)>"+bt2L+"&&(chtrg2r.ped-chtrg2r.min)<"+bt2R+" && track.trackis==1&&track.vz>=0.9998";
	//exbeamdata->Draw(DATA,CONDITIONwoCutsOnEn,"COLZ");
	exbeamdata->Draw(DATA,"","COLZ");

	DATA = "track.ys:track.xs>>hxyGEMswoXYcuts";
	TString CONDITIONWOGEMCUT = "(chtrg1r.ped-chtrg1r.min)>"+bt1L+"&&(chtrg1r.ped-chtrg1r.min)<"+bt1R+"&&(chtrg2r.ped-chtrg2r.min)>"+bt2L+"&&(chtrg2r.ped-chtrg2r.min)<"+bt2R+" &&lecroy2249.ch0>"+bCEnL+"&&lecroy2249.ch0<"+bCEnR+" && track.trackis==1 && track.vz>=0.9998";
	exbeamdata->Draw(DATA,CONDITIONWOGEMCUT,"colz");

	DATA = "track.ys:track.xs>>hxyGEMsXYWLS";
        TString CONDITIONGEMCUTWLS = "(chtrg1r.ped-chtrg1r.min)>"+bt1L+"&&(chtrg1r.ped-chtrg1r.min)<"+bt1R+"&&(chtrg2r.ped-chtrg2r.min)>"+bt2L+"&&(chtrg2r.ped-chtrg2r.min)<"+bt2R+"&& lecroy2249.ch0>"+bCEnL+"&&lecroy2249.ch0<"+bCEnR+" && track.trackis==1 && track.vz>=0.9998 && ch1.amp>1550";
        exbeamdata->Draw(DATA,CONDITIONGEMCUTWLS,"colz");

	//------Set selections-----------------------------
	CONDITION = "track.ys>"+by1+"   &&track.ys<"+by2+"  &&track.xs>"+bx1+"   &&track.xs<"+bx2+"  &&(chtrg1r.ped-chtrg1r.min)>"+bt1L+"&&(chtrg1r.ped-chtrg1r.min)<"+bt1R+"&&(chtrg2r.ped-chtrg2r.min)>"+bt2L+"&&(chtrg2r.ped-chtrg2r.min)<"+bt2R+"&&lecroy2249.ch0>"+bCEnL+"&&lecroy2249.ch0<"+bCEnR+" && track.trackis==1 && track.vz>=0.9998";
	
	//!!!!w/o GEM Cuts
	//CONDITION = "(chtrg1r.ped-chtrg1r.min)>"+bt1L+"&&(chtrg1r.ped-chtrg1r.min)<"+bt1R+"&&(chtrg2r.ped-chtrg2r.min)>"+bt2L+"&&(chtrg2r.ped-chtrg2r.min)<"+bt2R+"&&lecroy2249.ch0>"+bCEnL+"&&lecroy2249.ch0<"+bCEnR+" && track.trackis==1 && track.vz>=0.9998";
	
	DATA = "lecroy2249.ch0>>hcalCuts";
	hcalCuts->SetLineColor(kRed);
	exbeamdata->Draw(DATA,CONDITION,"same");

	DATA = "track.ys:track.xs>>hxyGEMs0";
	exbeamdata->Draw(DATA,CONDITION,"colz");
	
	DATA = "chtrg1r.ped-chtrg1r.min>>hchtrg1";
	exbeamdata->Draw(DATA,"","goff");

	DATA = "chtrg2r.ped-chtrg2r.min>>hchtrg2";
	hchtrg2->SetLineColor(7);
	exbeamdata->Draw(DATA,"","goff");

	//DATA="(ch1r.ped-ch1r.min):ch1.ti>>htimeAmp";
	DATA="ch1.amp:ch1.ti>>htimeAmp";
        exbeamdata -> Draw(DATA, CONDITION,"");

        //**********************************************
	//channel1
	//T=15
	//float A1phe = 13.527400; //U53.5V
	//float A1phe = 15.217; //U54V
	//float A1phe = 20.6155; //U55V
	//float A1phe = 27.7749; //U56V
	//T=45
	//float A1phe = 31.134; //U56V
	//float A1phe = 35.4146; //U57V
	//float A1phe = 50.0208; //U58V
	//
        float A1phe_int = 133.353;  //not use integral
      
	//channel2
	//T=15
	float A1phe = 13.495150; //U53.5V
	//float A1phe = 15.8427; //U54V
	//float A1phe = 21.0927; //U55V
	//float A1phe = 26.1076; //U56V
	//T=45
	//float A1phe = 22.0807; //42.322; //U56V    
	//float A1phe = 29.8345; // 47.1637; //U57V     
	//float A1phe = 37.6456; //60.465; //U58V	 
        //**********************************************
        Double_t binPed1 = 550.;
        Double_t binPed2 = 650.;

        Double_t binSig1 = 660.;
        Double_t binSig2 = 760.;

	///////
	Double_t binStarFit = 0.;

	vector<float> sigA; 	
	vector<float> sigIntA; 	
	vector<float> GEMX; 	
	vector<float> GEMY; 	
 
        for(Int_t i=0; i<exbeamdata->GetEntries(); i++)
        //for(Int_t i=0; i<50; i++)
	{
		binStarFit = binPed1;

		DATA="ch1.amp:ch1.ti >> amp0";
        	exbeamdata->Draw(DATA, CONDITION,"",1,i);
	
		if( exbeamdata->Draw(DATA, CONDITION,"",1,i)<=0 ) continue;

		auto h = new TH1F("h","osc",1023,0,1023);
		for (int j = 0; j < amp0->GetNbinsX(); j++ ) 
		{
            	    h->SetBinContent(j,amp0->GetBinContent(j)) ;
		    //cout<<j<<"\t"<<amp0->GetBinContent(j)<<endl;
        	}

		auto s = new TSpectrum();
        	auto hb = s->Background(h,40,"same");

		TProfile* amp = new TProfile("amp", "a(t) for Cnt", 1024, 0. , 1024);
		//cout<<h->GetNbinsX()<<endl;
	        for(int k=0; k<=h->GetNbinsX(); k++)
        	{
                   amp->Fill( k, h->GetBinContent(k) - hb->GetBinContent(k) );
		   //cout<<k<<"\t"<<h->GetBinContent(k)<<"\t"<<hb->GetBinContent(k)<<endl;
                }
		
	        //ped area
		amp->GetXaxis()->SetRangeUser(binPed1, binPed2);

		float meanPed = amp->GetMean(2);
	        float maxped1 = amp->GetMaximum();
        	float maxbinped1 = amp->GetMaximumBin();
                float maxped2 = amp->GetBinContent(maxbinped1+1);
         	float maxped3 = amp->GetBinContent(maxbinped1+2);
		float ped = (maxped1+maxped2+maxped3)/3.;

		float intPed = 0.;
        	for(int h=binPed1; h<=binPed2; h++)
                {
                   Double_t binContPed = amp->GetBinContent(h);
                   intPed += binContPed;
		}

		//signal area 
	        amp->GetXaxis()->SetRangeUser(binSig1, binSig2);

		float max1 = amp->GetMaximum();
	        float maxbin1 = amp->GetMaximumBin();
        	float max2 = amp->GetBinContent(maxbin1+1);
                float max3 = amp->GetBinContent(maxbin1+2);
		float sig = (max1+max2+max3)/3.;

		float intSig = 0.;
	        for(int n=binSig1; n<=binSig2; n++)
        	{
                   Double_t binContSig = amp->GetBinContent(n);
		   intSig += binContSig;
		}

		hAmpPed->Fill(ped);
		hAmpSig->Fill(sig);
		sigA.push_back(sig);
		 
		hAmpPedInt->Fill(intPed);
 	        hAmpSigInt->Fill(intSig);
		sigIntA.push_back(intSig);

                amp->GetXaxis()->SetRangeUser(0., 1024.);
	
		DATA = "track.ys:track.xs>>hxyGEMs";
		exbeamdata->Draw(DATA, CONDITION,"",1,i);
		//cout<<i<<"\t"<<hxyGEMs->GetMean(1)<<"\t"<<hxyGEMs->GetMean(2)<<endl;
		GEMX.push_back(hxyGEMs->GetMean(1));
		GEMY.push_back(hxyGEMs->GetMean(2));

		if( waveformdraw==1 && i<5000 )
		{
		  TString ii;
        	  ii.Form("%d", i);
	
		  amp0->Write("amp0_"+ii);
		  amp->Write("amp_"+ii);
		  hxyGEMs->Write("hxyGEMs_"+ii);
		}

	        for(int j=0; j<1024; j++)
        	{
                   Double_t binAmp2 = amp->GetBinContent(j);
                   htimeAmpSelected->Fill(j,binAmp2);
                }

		amp->Delete();
		h->Delete();
			
	}       	 
	
       	//fill hists with times distributions
	//-------------------------------
	if(TimeDist== 1)
	{
	DATA = "(Ttrig2-Ttrig1)>>hdtTrg";
	exbeamdata->Draw(DATA,CONDITION,"COLZ");

	TF1*G2 = new TF1("G2","[0]*exp(-(x-[1])^2/(2*[2]^2))+[3]*exp(-(x-[4])^2/(2*[5]^2))",9.5,11);
	G2->SetParameters(hdtTrg->GetMaximum(),hdtTrg->GetMean(),hdtTrg->GetRMS(),hdtTrg->GetMaximum()/5,hdtTrg->GetMean(),hdtTrg->GetRMS()*4);
	G2->SetParLimits(0,0,hdtTrg->GetMaximum()+10);
	G2->SetParLimits(1,hdtTrg->GetMean()-1*hdtTrg->GetRMS(),hdtTrg->GetMean()+1*hdtTrg->GetRMS());
	G2->SetParLimits(2,0,0.5);
	G2->SetParLimits(3,1,hdtTrg->GetMaximum()/1.05);
	G2->SetParLimits(4,hdtTrg->GetMean()-0.1*hdtTrg->GetRMS(),hdtTrg->GetMean()+0.1*hdtTrg->GetRMS());
	G2->SetParLimits(5,hdtTrg->GetRMS(),0.1);
	G2->SetParNames("Const_{#Delta T}","Mean_{#Delta T}","#sigma_{#Delta T}","Const_{dist}","Mean_{dist}","#sigma_{dist}");
	hdtTrg->Fit(G2,"E1","same",-10.0,10.0);
	hdtTrg->GetXaxis()->SetTitle("#Delta T, ns");
	hdtTrg->GetYaxis()->SetTitle("N");
	gStyle->SetOptFit(1111);

 	DATA = "(Ttrig1 - Tmcp)>>hdtCntTrg1";
	exbeamdata->Draw(DATA,CONDITION,"COLZ");

	TF1*G2_1 = new TF1("G2","[0]*exp(-(x-[1])^2/(2*[2]^2))+[3]*exp(-(x-[4])^2/(2*[5]^2))",9.5,11);
	G2_1->SetParameters(hdtCntTrg1->GetMaximum(),hdtCntTrg1->GetMean(),hdtCntTrg1->GetRMS(),hdtCntTrg1->GetMaximum()/5,hdtCntTrg1->GetMean(),hdtCntTrg1->GetRMS()*4);
	G2_1->SetParLimits(0,0,hdtCntTrg1->GetMaximum()+10);
	G2_1->SetParLimits(1,hdtCntTrg1->GetMean()-1*hdtCntTrg1->GetRMS(),hdtCntTrg1->GetMean()+1*hdtCntTrg1->GetRMS());
	G2_1->SetParLimits(2,0,0.5);
	G2_1->SetParLimits(3,1,hdtCntTrg1->GetMaximum()/1.05);
	G2_1->SetParLimits(4,hdtCntTrg1->GetMean()-0.1*hdtCntTrg1->GetRMS(),hdtCntTrg1->GetMean()+0.1*hdtCntTrg1->GetRMS());
	G2_1->SetParLimits(5,hdtCntTrg1->GetRMS(),0.1);
	G2_1->SetParNames("Const_{#Delta T}","Mean_{#Delta T}","#sigma_{#Delta T}","Const_{dist}","Mean_{dist}","#sigma_{dist}");
	hdtCntTrg1->Fit(G2_1,"E1","same",-10.0,10.0);
	hdtCntTrg1->GetXaxis()->SetTitle("#Delta T, ns");
	hdtCntTrg1->GetYaxis()->SetTitle("N");
	gStyle->SetOptFit(1111);      
 		
	DATA = "(Tmcp-Ttrig2)>>hdtCntTrg2";
	exbeamdata->Draw(DATA,CONDITION,"COLZ");

	TF1*G2_2 = new TF1("G2","[0]*exp(-(x-[1])^2/(2*[2]^2))+[3]*exp(-(x-[4])^2/(2*[5]^2))",9.5,11);
	G2_2->SetParameters(hdtCntTrg2->GetMaximum(),hdtCntTrg2->GetMean(),hdtCntTrg2->GetRMS(),hdtCntTrg2->GetMaximum()/5,hdtCntTrg2->GetMean(),hdtCntTrg2->GetRMS()*4);
	G2_2->SetParLimits(0,0,hdtCntTrg2->GetMaximum()+10);
	G2_2->SetParLimits(1,hdtCntTrg2->GetMean()-1*hdtCntTrg2->GetRMS(),hdtCntTrg2->GetMean()+1*hdtCntTrg2->GetRMS());
	G2_2->SetParLimits(2,0,0.5);
	G2_2->SetParLimits(3,1,hdtCntTrg2->GetMaximum()/1.05);
	G2_2->SetParLimits(4,hdtCntTrg2->GetMean()-0.1*hdtCntTrg2->GetRMS(),hdtCntTrg2->GetMean()+0.1*hdtCntTrg2->GetRMS());
	G2_2->SetParLimits(5,hdtCntTrg2->GetRMS(),0.1);
	G2_2->SetParNames("Const_{#Delta T}","Mean_{#Delta T}","#sigma_{#Delta T}","Const_{dist}","Mean_{dist}","#sigma_{dist}");
	hdtCntTrg2->Fit(G2_2,"E1","same",-10.0,10.0);
	hdtCntTrg2->GetXaxis()->SetTitle("#Delta T, ns");
	hdtCntTrg2->GetYaxis()->SetTitle("N");
	gStyle->SetOptFit(1111);

	DATA = "abs(Tmcp-Ttrig1):(ch1r.ped-ch1r.min)>>hdtCntTrg1AmpCnt";
	TString CONDITIONdtCntTrg1AmpCnt = "track.ys>"+by1+"   &&track.ys<"+by2+"  &&track.xs>"+bx1+"   &&track.xs<"+bx2+"  &&(chtrg1r.ped-chtrg1r.min)>"+bt1L+"&&(chtrg1r.ped-chtrg1r.min)<"+bt1R+"&&(chtrg2r.ped-chtrg2r.min)>"+bt2L+"&&(chtrg2r.ped-chtrg2r.min)<"+bt2R+"&& (ch1r.ped-ch1r.min)>0 &&lecroy2249.ch0>"+bCEnL+"&&lecroy2249.ch0<"+bCEnR+" &&track.trackis==1&&track.vz>=0.9998 && abs(Tmcp-Ttrig1)>0.0";
	exbeamdata->Draw(DATA,CONDITIONdtCntTrg1AmpCnt,"");
		
	TF1* fun2=new TF1("logt2","[0]+[1]*log(x+[2])",0,500);
	fun2->SetParLimits(1,0.0,100);
	hdtCntTrg1AmpCnt->Fit(fun2,"E1","same",10,550);
	hdtCntTrg1AmpCnt->GetXaxis()->SetTitle("A, ch");
	hdtCntTrg1AmpCnt->GetYaxis()->SetTitle("#Delta T, ns");
	gStyle->SetOptFit(11111);
        
	DATA = "abs(Tmcp-Ttrig2):(ch1r.ped-ch1r.min)>>hdtCntTrg2AmpCnt";
	TString CONDITIONdtCntTrg2AmpCnt = "track.ys>"+by1+"   &&track.ys<"+by2+"  &&track.xs>"+bx1+"   &&track.xs<"+bx2+"  &&(chtrg1r.ped-chtrg1r.min)>"+bt1L+"&&(chtrg1r.ped-chtrg1r.min)<"+bt1R+"&&(chtrg2r.ped-chtrg2r.min)>"+bt2L+"&&(chtrg2r.ped-chtrg2r.min)<"+bt2R+"&& (ch1r.ped-ch1r.min)>0 &&lecroy2249.ch0>"+bCEnL+"&&lecroy2249.ch0<"+bCEnR+" &&track.trackis==1&&track.vz>=0.999 && abs(Tmcp-Ttrig2)>0.0";
	exbeamdata->Draw(DATA,CONDITIONdtCntTrg2AmpCnt,"");
		
	TF1* fun=new TF1("logt","[0]+[1]*log(x+[2])",0,500);
	fun->SetParLimits(1,0.0,0.025);
	fun->SetParLimits(2,0,300);
	hdtCntTrg2AmpCnt->Fit(fun,"E1","same",10,550);
	hdtCntTrg2AmpCnt->GetXaxis()->SetTitle("A, ch");
	hdtCntTrg2AmpCnt->GetYaxis()->SetTitle("#Delta T, ns");
	gStyle->SetOptFit(11111);
	}
        //-----------------------------------------------------------

        //Npe calculation
        //------------------------	
	cout<<"1 method - Amplitudes:"<<endl;	
	
	//Float_t border1 = 13.5;
	hAmpPed->SetAxisRange(-16., 16.,"X");
	TF1* fitG = new TF1("fitG","gaus",-10.,50.);
        hAmpPed->Fit("fitG","LQ","",-10.,9.00);
        Double_t g1 = fitG->GetParameter(1);
        Double_t g2 = fitG->GetParameter(2);
        float border1 = g1 + 2.0*g2;
	cout<<"ApedMean ="<<hAmpPed->GetMean() <<" +/- "<< hAmpPed->GetRMS()/sqrt(hAmpPed->GetEntries())<<"\t"<<"ApedRMS="<<hAmpPed->GetRMS()<<endl;
	cout<<"From fit ApedMean ="<< g1 <<" +/- "<<"\t"<<"ApedRMS="<< g2 <<endl;
	cout<<"Board for pedestal in ch1 ampl spectr is "<<border1<<endl;

	cout<<"Ch_Ampl. = " << hAmpSig->GetMean() <<" +/- "<<hAmpSig->GetRMS()/sqrt(hAmpSig->GetEntries())<<endl;
        cout<<"Ch_Ampl. - Aped = " << hAmpSig->GetMean() - hAmpSig->GetMean() <<endl;
	cout<<"A1phe="<<A1phe<<endl;
	//cout<<"Ch_Ampl. - Aped = " << hAmpSig->GetMean() - g1 <<endl;
	cout<<"Aped Entries is "<<hAmpPed->GetEntries()<<"\t"<<"Ch_Ampl. Entries is "<<hAmpSig->GetEntries()<<endl;

	float pedVal = hAmpPed->GetMean();

	float resNpe = 0;
	float errresNpe = 0;
	//resNpe = ( hAmpSig->GetMean() - hAmpPed->GetMean() )/A1phe;     
	//errresNpe = ( hAmpSig->GetRMS()/sqrt(hAmpSig->GetEntries()) + hAmpPed->GetRMS()/sqrt(hAmpPed->GetEntries()) )/A1phe;
	resNpe = ( hAmpSig->GetMean() - g1 )/A1phe;     
	int firstbin = hAmpSig->GetXaxis()->FindBin(-100);
	int binBorderPed = hAmpSig->GetXaxis()->FindBin(g1+3*g2);
	float nevPed = hAmpSig->Integral(firstbin, binBorderPed);
	errresNpe = ( hAmpSig->GetRMS()/sqrt(hAmpSig->GetEntries()-nevPed) + g2/sqrt(nevPed) )/A1phe;
	
	cout<<"Npe="<<resNpe<<"+-"<<errresNpe<<endl;
	
	hAmpPed->SetAxisRange(-100., 1900.,"X");
	hAmpSig->SetAxisRange(-100., 1900.,"X");

	///////////////////////////////////////////////////////////////
	float resNpeNoise = 0;
	resNpeNoise = ( hAmpPed->GetMean() - g1 )/A1phe;     
	cout<<"NpeNoise(norm)="<<resNpeNoise<<endl;

	int NbinBorder1 = hAmpPed->GetXaxis()->FindBin(border1);
        int Nfirstbin = hAmpPed->GetXaxis()->FindBin(-100);
        int Nlastbin = hAmpPed->GetXaxis()->FindBin(1900);

        float Nintegral_before = hAmpPed->Integral(Nfirstbin,NbinBorder1);
        float Nintegral_after = hAmpPed->Integral(NbinBorder1,Nlastbin);
        float Nintegral_all = hAmpPed->Integral(Nfirstbin,Nlastbin);

        float Neffitiency = (Nintegral_after/Nintegral_all);
        float NeffitiencyErr = sqrt(Nintegral_after*(Nintegral_all-Nintegral_after)/pow(Nintegral_after,3));

        float NpeNoise = -log(1-Neffitiency);
        float NpeErrNoise = NeffitiencyErr/(1-Neffitiency);
	cout<<"NpeNoise(poisson)="<<NpeNoise<<"+-"<<NpeErrNoise<<endl;
        //////////////////////////////////////////////////////////// 


	TH3* hGEMNphe = new TH3D("hGEMNphe", "GEM-XY from Nphe", 120, 0, 120, 80, 0, 80, 100, 0.0, 100.0);
	
	TString outF = outDir + list + ".dat";
    	//ofstream fileOut;	
    	//fileOut.open(outF,ios_base::app | ios_base::ate);
    	//fileOut.open(outF,ios_base::out | ios_base::ate); //!!!

	for(int n=0; n<sigA.size(); n++)
        {
           //cout<< n <<"\t"<< GEMX[n] <<"\t"<< GEMY[n] <<"\t"<< (sigA[n] - g1)/A1phe <<endl;		
	   //hnpeAmpSig->Fill((sigA[n]-pedVal)/A1phe);	
	   hnpeAmpSig->Fill((sigA[n] - g1)/A1phe);	
	   hGEMNphe->Fill( GEMX[n], GEMY[n], (sigA[n] - g1)/A1phe );

	   //fileOut << GEMX[n] <<"\t"<< GEMY[n] <<"\t"<< (sigA[n] - g1)/A1phe << "\n";
        }
	//fileOut.close();
	
	hnpeAmpSig->Draw();
	hnpeAmpSig->SetOption("HIST");
	
	hGEMNphe->Draw();	
	hGEMNphe->SetOption("COLZ");	

	auto profhGEMNphe = hGEMNphe->Project3DProfile("yx");
	profhGEMNphe->Draw("colz");	
	
	//====================================================================
	cout<<"2 method - Integration:"<<endl;	
	
	float resNpe_int = ( hAmpSigInt->GetMean()-hAmpPedInt->GetMean() )/A1phe_int;
	float errresNpe_int = ( hAmpSigInt->GetRMS()/sqrt(hAmpSigInt->GetEntries()) + hAmpPedInt->GetRMS()/sqrt(hAmpPedInt->GetEntries()) )/A1phe_int;

	float pedValInt = hAmpPedInt->GetMean();

	for(int n2=0; n2<sigIntA.size(); n2++)
        {
	   hnpeIntSig->Fill((sigIntA[n2]-pedValInt)/A1phe_int);	
        }
	hnpeIntSig->Draw();
	hnpeIntSig->SetOption("HIST");

	cout<<"Int_Ped ="<<hAmpPedInt->GetMean()<<endl;
	cout<<"Int_Sig - Int_Ped = " << hAmpSigInt->GetMean() - hAmpPedInt->GetMean() <<endl;
	cout<<"Npe_int="<<resNpe_int<<"+-"<<errresNpe_int<<endl;

        //-----------------------------------------------------------
	//write results in root-file
        hxyGEMs0->Write();
        hxyGEMswoXYcuts->Write();
	hxyGEMsXYWLS->Write();
        hchtrg1->Write();
        hchtrg2->Write();
	hAmpSig->SetLineColor(kRed);
        hAmpSig->Write();
	hAmpPed->SetLineColor(kBlue);
        hAmpPed->Write();
	hAmpPedInt->SetLineColor(kBlue);
        hAmpPedInt->Write();
	hAmpSigInt->SetLineColor(kRed);
        hAmpSigInt->Write();
	hnpeAmpSig->Write();
	hnpeIntSig->Write();
	htimeAmp->Write();
	htimeAmpSelected->Write();
	hcal->Write();
	hcalCuts->Write();
	hGEMNphe->Write();
	profhGEMNphe->Write();
	if(TimeDist== 1)
        {
	  hdtCntTrg2->Write();
          hdtTrg->Write();
          hdtCntTrg1->Write();
          hdtCntTrg1AmpCnt->Write();
          hdtCntTrg2AmpCnt->Write();
	}
        f->Close();
}
