// Make milliqan validation plots. D. Stuart, Oct. 2017.
// adapted F. Setti, Jul. 2019.
#include "/home/users/fsetti/milliQan/Headers/BmuonAnalysis.cc"

void makeProbGraph( ) 
{
	TCanvas *c = new TCanvas("c","A Simple Graph Example",800,600);
	float x[7] = { 0.005, 0.01, 0.02 , 0.03, 0.05, 0.07, 0.1 };
        float yOne[7] = { 0.105174, 0.352006, 0.79802, 0.965229, 0.999815, 1., 1. };
	float yNo[7]  = {0.894826, 0.647994,0.20198,0.0347715,0.000184963,0};
	float yAll[7] = {8.04699e-5,0.00187541,0.024285,0.141222,0.687043,0.956035,0.999253};
	
	TGraph *grOne, *grNo, *grAll;
	grOne = new TGraph (7,x,yOne);
	grNo  = new TGraph (7,x,yNo );
	grAll = new TGraph (7,x,yAll);

	grOne->SetMarkerStyle(17);
	grNo->SetMarkerStyle(17);
	grAll->SetMarkerStyle(17);

	grOne->SetLineWidth(2);
	grNo->SetLineWidth(2);
	grAll->SetLineWidth(2);
	
	grOne->SetLineColor(kGreen);
	grNo->SetLineColor(kBlue);
	grAll->SetLineColor(kOrange);

	grOne->SetMaximum(3.);
	grOne->SetMinimum(0.01);

	c->SetLogy();
	c->SetGridy();
	c->SetLogx();
	c->SetGridx();
	grOne->Draw("AC*");
	grNo->Draw("C*");
	grAll->Draw("C*");

	c->SaveAs("/home/users/fsetti/public_html/milliQan/SimDataComparison/14Oct2019/probGraph.png");
}
