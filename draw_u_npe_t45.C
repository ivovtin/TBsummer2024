void draw_u_npe_t45()
{
   double U[3] = {56, 57, 58};
   //point2
   //double npe[3] = {5.72, 6.57, 6.12};
   //double err_npe[3] = {0.10, 0.10, 0.10};
   //
   //point4
   double npe[3] = {5.00, 5.78, 5.40};
   double err_npe[3] = {0.10, 0.09, 0.10};
   //
   TCanvas *c1 = new TCanvas();
   c1->cd();	
   auto g = new TGraphErrors(3,U,npe,0,err_npe);

   //g->SetTitle("Point 2; U, V; N_{ph.e.}");
   g->SetTitle("Point 4; U, V; N_{ph.e.}");

   g->SetMarkerStyle(21);
   //g->SetMarkerColor(4);
   g->Draw();
   //c1->SaveAs("results/point2_U-Nphe.png");
   c1->SaveAs("results/point4_U-Nphe.png");
} 
