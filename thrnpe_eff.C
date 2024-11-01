void thrnpe_eff()
{
	//ch1
	TString infile = "ch1_point2_U53V_T15_res.root";	
	//TString infile = "ch1_point9_U53V_T15_res.root";	
	//TString infile = "ch1_point2_U54V_T15_res.root";	
	//TString infile = "ch1_point2_U55V_T15_res.root";	
	//TString infile = "ch1_point2_U56V_T15_res.root";	
	//TString infile = "ch1_point2_U57V_T45_res.root";	
	//TString infile = "ch1_point2_U56V_T45_res.root";	
	//TString infile = "ch1_point2_U58V_T45_res.root";	
	//TString infile = "ch1_point9_U57V_T45_res.root";	
	
	//ch2
	//TString infile = "ch2_point2_U53V_T15_res.root";	
	//TString infile = "ch2_point9_U53V_T15_res.root";	
	//TString infile = "ch2_point2_U57V_T45_res.root";	
	//TString infile = "ch2_point9_U57V_T45_res.root";	
	
	TString outDir="results/";
	TString outfile = "underther_eff_" + infile;

	TString prefix = "good_res_x61-68_y27-37/";

        //**********************************************
        //channel1
        //T=15
        float A1phe = 13.527400; //U53.5V
        //float A1phe = 15.217; //U54V
        //float A1phe = 20.6155; //U55V
        //float A1phe = 27.7749; //U56V
        //T=45
        //float A1phe = 31.134; //U56V
        //float A1phe = 35.4146; //U57V
        //float A1phe = 50.0208; //U58V

        //channel2
        //T=15
        //float A1phe = 13.495150; //U53.5V
        //float A1phe = 15.8427; //U54V
        //float A1phe = 21.0927; //U55V
        //float A1phe = 26.1076; //U56V
        //T=45
        //float A1phe = 22.0807; //42.322; //U56V     
        //float A1phe = 29.8345; //47.1637; //U57V    
        //float A1phe = 37.6456; //60.465; //U58V      
        //**********************************************

	TFile *f=new TFile(outDir + infile);
	//TFile *f=new TFile(outDir + prefix + infile);

	auto fout = TFile::Open(outDir + outfile,"RECREATE");
        TH1F *hped=(TH1F*)f->Get("hAmpPed");

	TF1* fitG = new TF1("fitG","gaus",-10.,50.);
        hped->Fit("fitG","LQ","",-10.,10.00);
        Double_t g1 = fitG->GetParameter(1);

        TH1F *hsig=(TH1F*)f->Get("hAmpSig");
	//hped->Draw();	
	
	hped->SetAxisRange(-10., 100.,"X");

	int binBorder1, firstbin, lastbin;
	float integral_before, integral_after, integral_all;
	float underthr_eff, underthr_effErr, underthr_Npe, underthr_NpeErr;
	
	int binBorder2, firstbin2, lastbin2;
	float integral_before2, integral_after2, integral_all2;
	float eff, effErr, Npe, NpeErr;
	float Ksigma, errKsigma, errKsigma1, errKsigma2;

	float border1;
	
	vector< float > vec_ch, vec_npe, vec_underthr_eff, vec_underthr_effErr, vec_eff, vec_effErr, vec_Ksigma, vec_errKsigma;

	for(int i=0; i<=110; i=i+2 ){    //p2, 53.5 t=15
	//for(int i=9; i<=54; i++ ){    //p9, 53.5 t=15
	//for(int i=9; i<=110; i=i+2 ){    //p2, 54 t=15
	//for(int i=9; i<=95; i++ ){    //p2, 55 t=15
	//for(int i=9; i<=95; i++ ){    //p2, 56 t=15
	//for(int i=31; i<=204; i+=5 ){    //p2, 57 t=45
	//for(int i=31; i<=204; i+=5 ){    //p9, 57 t=45
	//for(int i=31; i<=130; i+=5 ){    //p2, 56 t=45
	//for(int i=35; i<=270; i+=5 ){    //p2, 58 t=45
	//ch2
	//for(int i=19; i<=49; i++ ){    //p2, 53.5 t=15
	//for(int i=19; i<=47; i++ ){    //p9, 53.5 t=15
	//for(int i=52; i<=180; i+=5 ){    //p2, 57 t=45
	//for(int i=31; i<=203; i+=5 ){    //p9, 57 t=45
		cout<<"===================================================="<<endl;
		border1 = i;
		cout<<"border1 = "<<border1<<endl;
        	binBorder1 = hped->GetXaxis()->FindBin(border1);
        	firstbin = hped->GetXaxis()->FindBin(-100);
        	lastbin = hped->GetXaxis()->FindBin(1900);

        	integral_before = hped->Integral(firstbin,binBorder1);
        	integral_after = hped->Integral(binBorder1,lastbin);
        	integral_all = hped->Integral(firstbin,lastbin);
        	//cout<<"integral_after:"<<integral_after<<endl;
        	//cout<<"integral_all:"<<integral_all<<endl;

    	        if( integral_after<=1 )
                {
                        integral_after = 2;
                        integral_all = integral_all - 2;
                }

        	underthr_eff = (integral_after/integral_all);
		underthr_effErr = sqrt(integral_after*(integral_all-integral_after)/pow(integral_all,3));
        	underthr_Npe = -log(1-underthr_eff);
        	underthr_NpeErr = underthr_effErr/(1-underthr_eff);

        	binBorder2 = hsig->GetXaxis()->FindBin(border1);
        	firstbin2 = hsig->GetXaxis()->FindBin(-100);
        	lastbin2 = hsig->GetXaxis()->FindBin(1900);

        	integral_before2 = hsig->Integral(firstbin2,binBorder2);
        	integral_after2 = hsig->Integral(binBorder2,lastbin2);
        	integral_all2 = hsig->Integral(firstbin2,lastbin2);
        	//cout<<"integral_after2:"<<integral_after2<<endl;
        	//cout<<"integral_all2:"<<integral_all2<<endl;

		if( integral_after2>=(integral_all2-1) )
                {
                        integral_after2 = integral_after2 - 2;
                }

        	eff = (integral_after2/integral_all2);
		effErr = sqrt(integral_after2*(integral_all2-integral_after2)/pow(integral_all2,3));
        	Npe = -log(1-eff);
       	 	NpeErr = effErr/(1-eff);

		vec_ch.push_back(border1);
		//vec_npe.push_back(border1/A1phe);
		//vec_npe.push_back((border1 - hped->GetMean())/A1phe);
		vec_npe.push_back((border1 - g1)/A1phe);
		vec_underthr_eff.push_back(underthr_eff);
		vec_underthr_effErr.push_back(underthr_NpeErr);
		//vec_eff.push_back(eff);
		vec_eff.push_back(1.-eff);
		vec_effErr.push_back(effErr);

		Ksigma = abs(sqrt(2.)*(TMath::ErfInverse(1-2*underthr_eff)+TMath::ErfInverse(1-2*(1-eff))));
		//cout<<"underthr_eff-underthr_effErr = "<<underthr_eff-underthr_effErr<<endl;
		//cout<<"eff-effErr = "<<eff-effErr<<endl;
		//cout<<"underthr_eff+underthr_effErr = "<<underthr_eff+underthr_effErr<<endl;
		//cout<<"eff+effErr = "<<eff+effErr<<endl;
		errKsigma1 = Ksigma - abs(sqrt(2.)*(TMath::ErfInverse(1-2*(underthr_eff-underthr_effErr))+TMath::ErfInverse(1-2*(1-(eff+effErr)))));
		errKsigma2 = Ksigma - abs(sqrt(2.)*(TMath::ErfInverse(1-2*(underthr_eff+underthr_effErr))+TMath::ErfInverse(1-2*(1-(eff-effErr)))));
		cout<<"Ksigma = "<<Ksigma<<endl;
		cout<<"errKsigma1 = "<<errKsigma1<<endl;
		cout<<"errKsigma2 = "<<errKsigma2<<endl;
		errKsigma = abs(errKsigma1);
		if( abs(errKsigma2)>abs(errKsigma) ) errKsigma = abs(errKsigma2);
		vec_Ksigma.push_back(Ksigma);
		vec_errKsigma.push_back(errKsigma);

        	cout<<"*************************************************"<<endl;
        	cout<<"underthr_eff = " << underthr_eff * 100 <<" +/- "<< underthr_effErr*100 <<"%"<<endl;
        	cout<<"underthr_Nph.e. = " << underthr_Npe <<" +/- "<< underthr_NpeErr <<endl;
        	cout<<"*************************************************"<<endl;
		cout<<"eff = " << eff * 100 <<" +/- "<< effErr*100 <<"%"<<endl;
        	cout<<"Nph.e. = " << Npe <<" +/- "<<NpeErr<<endl;
        	cout<<"*************************************************"<<endl;
	}

	TGraphErrors* gr1 = new TGraphErrors(vec_ch.size(), &vec_ch[0], &vec_underthr_eff[0], 0, &vec_underthr_effErr[0]);
	//TGraphErrors* gr1 = new TGraphErrors(vec_ch.size(), &vec_ch[0], &vec_underthr_eff[0], 0, 0);
	gr1->SetMarkerStyle(21);
  	gr1->SetMarkerSize(0.8);
  	gr1->SetMarkerColor(4);
  	gr1->SetLineWidth(2);
  	gr1->SetLineColor(4);

	TGraphErrors* gr2 = new TGraphErrors(vec_ch.size(), &vec_ch[0], &vec_eff[0], 0, &vec_effErr[0]);
	//TGraphErrors* gr2 = new TGraphErrors(vec_ch.size(), &vec_ch[0], &vec_eff[0], 0, 0);
	gr2->SetMarkerStyle(24);
  	gr2->SetMarkerSize(0.8);
  	gr2->SetMarkerColor(2);
  	gr2->SetLineWidth(2);
  	gr2->SetLineColor(2);

	TMultiGraph *mg = new TMultiGraph();
   	mg->Add(gr1);
   	mg->Add(gr2);
	mg->SetTitle("; Threshold, Channels; Misidentification");
	mg->GetYaxis()->SetRangeUser(0, 1);

	TGraphErrors* gr3 = new TGraphErrors(vec_ch.size(), &vec_ch[0], &vec_Ksigma[0], 0, &vec_errKsigma[0]);
	gr3->SetMarkerStyle(24);
  	gr3->SetMarkerSize(0.8);
  	gr3->SetMarkerColor(6);
  	gr3->SetLineWidth(2);
  	gr3->SetLineColor(6);
	gr3->SetTitle("; Threshold, Channels; #sigma");
	gr3->GetYaxis()->SetRangeUser(0, 10);

	//***********************************************
	
	TGraphErrors* gr4 = new TGraphErrors(vec_npe.size(), &vec_npe[0], &vec_underthr_eff[0], 0, &vec_underthr_effErr[0]);
	gr4->SetMarkerStyle(21);
  	gr4->SetMarkerSize(0.8);
  	gr4->SetMarkerColor(4);
  	gr4->SetLineWidth(2);
  	gr4->SetLineColor(4);

	TGraphErrors* gr5 = new TGraphErrors(vec_npe.size(), &vec_npe[0], &vec_eff[0], 0, &vec_effErr[0]);
	gr5->SetMarkerStyle(24);
  	gr5->SetMarkerSize(0.8);
  	gr5->SetMarkerColor(2);
  	gr5->SetLineWidth(2);
  	gr5->SetLineColor(2);

	TMultiGraph *mg2 = new TMultiGraph();
   	mg2->Add(gr4);
   	mg2->Add(gr5);
	mg2->SetTitle("; Threshold, npe; Misidentification");
	mg2->GetYaxis()->SetRangeUser(0, 1);

	TGraphErrors* gr6 = new TGraphErrors(vec_npe.size(), &vec_npe[0], &vec_Ksigma[0], 0, &vec_errKsigma[0]);
	gr6->SetMarkerStyle(24);
  	gr6->SetMarkerSize(0.8);
  	gr6->SetMarkerColor(6);
  	gr6->SetLineWidth(2);
  	gr6->SetLineColor(6);
	gr6->SetTitle("; Threshold, npe; #sigma");
	gr6->GetYaxis()->SetRangeUser(0, 10);

	mg->Write("misident");
	mg2->Write("misident_npe");
	gr3->Write("sigma");
	gr6->Write("sigma_npe");
	hped->SetLineColor(4);
        hped->SetTitle(" ");
	hped->Write();
	hsig->SetLineColor(2);
        hsig->SetTitle(" ");
	hsig->Write();
	fout->Close();
}	

