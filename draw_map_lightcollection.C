void draw_map_lightcollection()
{
	TString inDir = "/home/ovtin/development/TB_SND_ASHIPH/TBsummer2024/results/";
	float x, y, npe;

	float dx = 0.;
	float dy = 0.;

	TString outfile = "ch1_lightcollection.root";

	auto f = TFile::Open(inDir + outfile,"RECREATE");

	TH3* hGEMNphe = new TH3D("hGEMNphe", "Nphe from GEM-XY", 200, 0, 200, 200, 0, 250, 100, 0.0, 100);
	//TH2D *hGEMNphe = new TH2D("hGEMNphe", "Nphe from GEM-XY", 50, 0, 200, 50, 0, 250);

	//TProfile2D *hGEMNphe  = new TProfile2D("hGEMNphe","Nphe from GEM-XY",200,0,200,200,0,250,0,100);

	//auto histx = new TH2F("histx", "", 200, 0, 200, 100, 0, 100);
	//auto histy = new TH2F("histy", "", 200, 0, 250, 100, 0, 100);

	auto profx = new TProfile("profx", "", 200, 0, 200, 0, 100);
	auto profy = new TProfile("profy", "", 200, 0, 250, 0, 100);


	for (string list : {"ch1_point1_U53V_T15.dat", "ch1_point2_U53V_T15.dat", "ch1_point3_U53V_T15.dat", "ch1_point4_U53V_T15.dat", "ch1_point5_U53V_T15.dat", 
			"ch1_point6_U53V_T15.dat", "ch1_point7_U53V_T15.dat", "ch1_point8_U53V_T15.dat", "ch1_point9_U53V_T15.dat", "ch1_point10_U53V_T15.dat",
		       "ch1_point11_U53V_T15.dat", "ch1_point12_U53V_T15.dat" })
	{

	    size_t point1 = list.find("point1");
	    size_t point2 = list.find("point2");
	    size_t point3 = list.find("point3");
	    size_t point4 = list.find("point4");
	    size_t point5 = list.find("point5");
	    size_t point6 = list.find("point6");
	    size_t point7 = list.find("point7");
	    size_t point8 = list.find("point8");
	    size_t point9 = list.find("point9");
	    size_t point10 = list.find("point10");
	    size_t point11 = list.find("point11");
	    size_t point12 = list.find("point12");

	    //cout<<point1<<endl;

	    if( point1==4 ){ dx = 0.; dy = 161.; }
	    if( point2==4 ){ dx = 25.; dy = 161.; }
	    if( point3==4 ){ dx = 47.; dy = 161.; }
	    if( point4==4 ){ dx = 67.; dy = 161.; }
	    if( point5==4 ){ dx = 0.; dy = 80.; }
	    if( point6==4 ){ dx = 25.; dy = 80.; }
	    if( point7==4 ){ dx = 47.; dy = 80.; }
	    if( point8==4 ){ dx = 67.; dy = 80.; }
	    if( point9==4 ){ dx = 0.; dy = 0.; }
	    if( point10==4 ){ dx = 25.; dy = 0.; }
	    if( point11==4 ){ dx = 47.; dy = 0.; }
	    if( point12==4 ){ dx = 67.; dy = 0.; }

	    ifstream fin(inDir + list);
	    string line;
	    if( !fin.is_open() )
	    {
	      cout << "File not open!!!\n";
 	    }
	    else
	    {
	       while( getline(fin, line) )
	       {
		  if( line.size()>=1 && line[0]!='#' )
		  {
		    std::istringstream i_str(line);
		    i_str >> x >> y >> npe ;
		    //cout<< x <<"\t"<< y <<"\t"<< npe <<endl;
		    hGEMNphe->Fill( x + dx, y + dy, npe );
		    if( (x + dx)<110. || (x + dx)>130. )
		    {
		    profx->Fill(x + dx, npe);
		    profy->Fill(y + dy, npe);
		    }
		  }
	       }
	       fin.close();
	    }

	}

	gStyle->SetOptStat(0);
	/*
	Int_t palette[8];
   	palette[0] = 1;
   	palette[1] = 4;
   	palette[2] = 7;
   	palette[3] = 10;
   	palette[4] = 15;
   	palette[5] = 20;
   	palette[6] = 25;
   	palette[7] = 50;
   	gStyle->SetPalette(8,palette);
	*/
	gStyle->SetPalette(1);
	//hGEMNphe->SetContour(10);
	hGEMNphe->SetContour(8);
	//hGEMNphe->GetZaxis()->SetRangeUser(0, 8);
	//hGEMNphe->GetZaxis()->SetRangeUser(0, 10);
	//hGEMNphe->GetZaxis()->SetRangeUser(10, 20);
	//hGEMNphe->GetZaxis()->SetRangeUser(20, 30);
	//hGEMNphe->GetZaxis()->SetRangeUser(20, 50);
	//hGEMNphe->GetZaxis()->SetRangeUser(30, 50);
	hGEMNphe->GetXaxis()->SetTitle("X, mm");
	hGEMNphe->GetYaxis()->SetTitle("Y, mm");
	hGEMNphe->Draw("COLZ");

	auto profhGEMNphe = hGEMNphe->Project3DProfile("yx");
   	profhGEMNphe->SetContour(12);
	profhGEMNphe->GetZaxis()->SetRangeUser(8, 20);
	profhGEMNphe->GetXaxis()->SetTitle("X, mm");
        profhGEMNphe->GetYaxis()->SetTitle("Y, mm");
   	profhGEMNphe->Draw("colz");

	profx->GetXaxis()->SetTitle("X, mm");	
	profx->GetYaxis()->SetTitle("N_{ph.e.}");
	//profx->Draw();
	profy->GetXaxis()->SetTitle("Y, mm");	
	profy->GetYaxis()->SetTitle("N_{ph.e.}");
	profy->Draw();

	hGEMNphe->Write();
	profhGEMNphe->Write();
	profx->Write();
	profy->Write();	
}
