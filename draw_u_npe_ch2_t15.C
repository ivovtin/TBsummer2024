void draw_u_npe_ch2_t15()
{
   double U[3] = {53.5, 54, 55};
   //point2
   double npe[3] = {12.71, 13.23, 13.30};
   double err_npe[3] = {0.19, 0.11, 0.12};
   //
   TCanvas *c1 = new TCanvas();
   c1->cd();	
   auto g = new TGraphErrors(3,U,npe,0,err_npe);

   g->SetTitle("Point 2; U, V; N_{ph.e.}");
   //g->SetTitle("Point 6; U, V; N_{ph.e.}");

   g->SetMarkerStyle(21);
   //g->SetMarkerColor(4);
   g->Draw();
   c1->SaveAs("results/point2_U-Nphe.png");
   //c1->SaveAs("results/point6_U-Nphe.png");
} 
