

int DefinePtBins(){
  float width[]={20-15, 26-20, 35-26, 46-35, 61-46, 80-61, 106-80, 141-106, 186-141, 247-186, 326-247, 432-326, 571-432, 756-571, 1000-756};
  float start[]={15, 20, 26, 35, 46, 61, 80, 106, 141, 186, 247, 326, 432, 571, 756};

  TGraph* gr = new TGraph(15,start,width);
  TCanvas * c1 = new TCanvas("c1","c1",800,600);
  gr->Draw("A*");
  c1->SaveAs("ptbin_width.png");

  return 0;
}
