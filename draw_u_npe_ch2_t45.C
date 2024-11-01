void draw_u_npe_ch2_t45()
{
   double U[3] = {56, 57, 58};
   //point2
   //double npe[3] = {12.4874, 11.97, 12.15};
   //double err_npe[3] = {0.10, 0.10, 0.10};
   //
   //point9
   double npe[3] = {7.61, 7.65, 7.73};
   double err_npe[3] = {0.10, 0.09, 0.09};
   //
   TCanvas *c1 = new TCanvas();
   c1->cd();	
   auto g = new TGraphErrors(3,U,npe,0,err_npe);

   //g->SetTitle("Point 2; U, V; N_{ph.e.}");
   g->SetTitle("Point 9; U, V; N_{ph.e.}");

   g->SetMarkerStyle(21);
   //g->SetMarkerColor(4);
   g->Draw();
   //c1->SaveAs("results/point2_U-Nphe.png");
   c1->SaveAs("results/point9_U-Nphe.png");
} 
