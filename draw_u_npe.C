void draw_u_npe()
{
   double U[4] = {53.5, 54, 55, 56};
   //point2
   double npe[4] = {7.14, 8.36, 9.16, 8.72};
   double err_npe[4] = {0.11, 0.08, 0.08, 0.08};
   //
   //point6
   //double npe[4] = {6.53, 7.12, 7.35, 7.37};
   //double err_npe[4] = {0.08, 0.10, 0.07, 0.07};
   //
   TCanvas *c1 = new TCanvas();
   c1->cd();	
   auto g = new TGraphErrors(4,U,npe,0,err_npe);

   g->SetTitle("Point 2; U, V; N_{ph.e.}");
   //g->SetTitle("Point 6; U, V; N_{ph.e.}");

   g->SetMarkerStyle(21);
   //g->SetMarkerColor(4);
   g->Draw();
   c1->SaveAs("results/point2_U-Nphe.png");
   //c1->SaveAs("results/point6_U-Nphe.png");
} 
