#include "tdrstyle.C"

using namespace std;

void Draw_single1Dhist(TH1F *hist, bool norm1, bool logy, const char *XaxisName, const char *Ytitle_unit, int digitafterdecimal, const char *legtext, const char *outname)
{
    gStyle->SetOptStat(0);

    hist->Sumw2();

    if (norm1 == true)
    {
        hist->Scale(1. / hist->Integral(-1, -1));
    }

    TCanvas *c = new TCanvas("c1", "", 500, 500);
    if (logy)
        c->SetLogy();
    c->cd();
    gPad->SetRightMargin(0.05);
    gPad->SetTopMargin(0.08);
    gPad->SetLeftMargin(0.13);
    //
    Float_t binwidth = hist->GetXaxis()->GetBinWidth(1);
    hist->GetXaxis()->SetTitle(XaxisName);
    if (norm1 == true)
    {
        if (digitafterdecimal == 1)
        {
            if (strcmp(Ytitle_unit, "") == 0)
                hist->GetYaxis()->SetTitle(Form("Normalized entries/(%.1f)", binwidth));
            else
                hist->GetYaxis()->SetTitle(Form("Normalized entries/(%.1f %s)", binwidth, Ytitle_unit));
        }
        else if (digitafterdecimal == 2)
        {
            if (strcmp(Ytitle_unit, "") == 0)
                hist->GetYaxis()->SetTitle(Form("Normalized entries/(%.2f)", binwidth));
            else
                hist->GetYaxis()->SetTitle(Form("Normalized entries/(%.2f %s)", binwidth, Ytitle_unit));
        }
        else if (digitafterdecimal == 3)
        {
            if (strcmp(Ytitle_unit, "") == 0)
                hist->GetYaxis()->SetTitle(Form("Normalized entries/(%.3f)", binwidth));
            else
                hist->GetYaxis()->SetTitle(Form("Normalized entries/(%.3f %s)", binwidth, Ytitle_unit));
        }
        else
        {
            cout << "To much digits! Use 1 digit after decimal point!" << endl;
            if (strcmp(Ytitle_unit, "") == 0)
                hist->GetYaxis()->SetTitle(Form("Normalized entries/(%.1f)", binwidth));
            else
                hist->GetYaxis()->SetTitle(Form("Normalized entries/(%.1f %s)", binwidth, Ytitle_unit));
        }
        
    }
    else
    {
        if (digitafterdecimal == 1)
        {
            if (strcmp(Ytitle_unit, "") == 0)
                hist->GetYaxis()->SetTitle(Form("Entries/(%.1f)", binwidth));
            else
                hist->GetYaxis()->SetTitle(Form("Entries/(%.1f %s)", binwidth, Ytitle_unit));
        }
        else if (digitafterdecimal == 2)
        {
            if (strcmp(Ytitle_unit, "") == 0)
                hist->GetYaxis()->SetTitle(Form("Entries/(%.2f)", binwidth));
            else
                hist->GetYaxis()->SetTitle(Form("Entries/(%.2f %s)", binwidth, Ytitle_unit));
        }
        else if (digitafterdecimal == 3)
        {
            if (strcmp(Ytitle_unit, "") == 0)
                hist->GetYaxis()->SetTitle(Form("Entries/(%.3f)", binwidth));
            else
                hist->GetYaxis()->SetTitle(Form("Entries/(%.3f %s)", binwidth, Ytitle_unit));
        }
        else
        {
            cout << "To much digits! Use 1 digit after decimal point!" << endl;
            if (strcmp(Ytitle_unit, "") == 0)
                hist->GetYaxis()->SetTitle(Form("Entries/(%.1f)", binwidth));
            else
                hist->GetYaxis()->SetTitle(Form("Entries/(%.1f %s)", binwidth, Ytitle_unit));
        }
        
    }

    Float_t rangescale = (logy) ? 10 : 1.5;
    if (logy == true)
        hist->GetYaxis()->SetRangeUser(1, (hist->GetMaximum()) * rangescale);
    else
        hist->GetYaxis()->SetRangeUser(0, (hist->GetMaximum()) * rangescale);

    hist->GetXaxis()->SetTickSize(0.03);
    hist->GetXaxis()->SetTitleSize(0.04);
    hist->GetXaxis()->SetLabelSize(0.03);
    hist->GetYaxis()->SetTickSize(0.03);
    hist->GetYaxis()->SetTitleSize(0.04);
    hist->GetYaxis()->SetLabelSize(0.03);
    hist->GetXaxis()->SetTitleOffset(1.3);
    hist->GetYaxis()->SetTitleOffset(1.3);
    hist->SetLineColor(1);
    hist->SetLineWidth(2);
    hist->Draw("hist");
    TLegend *l = new TLegend(0.5, 0.8, 0.89, 0.89);
    l->AddEntry(hist, legtext, "l");
    l->SetFillColor(0); //Set the background to be white
    l->SetLineColor(1);
    l->Draw("same");
    c->SaveAs(Form("%s.png", outname));
    c->SaveAs(Form("%s.pdf", outname));
}

void Gen_study(const char *infile)
{
    setTDRStyle();

    fstream fp, fp2;
    // fp.open("ParentInfo_PID13_hardproc_isprompt.csv", ios::out | ios::trunc);
    fp.open("ParentInfo_PID13_hardproc_isprompt_mZeq2mu.csv", ios::out | ios::trunc);
    fp << "ev,genmu1PID,genmu1MomPID,genmu1GMomPID,genmu2PID,genmu2MomPID,genmu2GMomPID" << endl;

    // fp2.open("ParentInfo_PID13.csv", ios::out | ios::trunc);
    // fp2 << "ev,genmuPID,genmuMomPID,genmuGMomPID" << endl;

    TH1F *hM_genMuPt = new TH1F("hM_genMuPt", "hM_genMuPt", 50, 0, 100);
    TH1F *hM_genMuEta = new TH1F("hM_genMuEta", "hM_genMuEta", 50, -5, 5);
    TH1F *hM_genMuPhi = new TH1F("hM_genMuPhi", "hM_genMuPhi", 80, -4, 4);
    TH1F *hM_genMmm = new TH1F("hM_genMmm", "hM_genMmm", 100, 0, 200);
    TH1F *hM_genMomMass = new TH1F("hM_genMomMass", "hM_genMomMass", 100, 0, 200);
    TH1F *hM_genMomMass_lowM = new TH1F("hM_genMomMass_lowM", "hM_genMomMass_lowM", 50, 0, 2);

    TFile *fin = new TFile(infile, "READ");
    TTree *tree = (TTree *)fin->Get("EventTree");

    Int_t nGen;
    vector<float> *genPt = 0;
    vector<float> *genEta = 0;
    vector<float> *genPhi = 0;
    vector<int> *genPID = 0;
    vector<int> *genMomPID = 0;
    vector<int> *genGMomPID = 0;
    vector<float> *genMomMass = 0;
    vector<UShort_t> *genStatusFlag = 0;
    Int_t nGenJet;
    vector<float> *genJetPt = 0;
    vector<float> *genJetEta = 0;
    vector<float> *genJetPhi = 0;
    vector<float> *genJetE = 0;

    tree->SetBranchAddress("nGen", &nGen);
    tree->SetBranchAddress("genPt", &genPt);
    tree->SetBranchAddress("genEta", &genEta);
    tree->SetBranchAddress("genPhi", &genPhi);
    tree->SetBranchAddress("genPID", &genPID);
    tree->SetBranchAddress("genMomPID", &genMomPID);
    tree->SetBranchAddress("genGMomPID", &genGMomPID);
    tree->SetBranchAddress("genMomMass", &genMomMass);
    tree->SetBranchAddress("genStatusFlag", &genStatusFlag);
    tree->SetBranchAddress("nGenJet", &nGenJet);
    tree->SetBranchAddress("genJetPt", &genJetPt);
    tree->SetBranchAddress("genJetEta", &genJetEta);
    tree->SetBranchAddress("genJetPhi", &genJetPhi);
    tree->SetBranchAddress("genJetE", &genJetE);

    for (Int_t ev = 0; ev < tree->GetEntries(); ++ev)
    {
        tree->GetEntry(ev);

        vector<int> vec_Mu;
        vec_Mu.clear();
        for (Int_t igen = 0; igen < nGen; ++igen)
        {
            if (abs(genPID->at(igen)) != 13)
                continue;

            // fp2 << ev << "," << genPID->at(igen) << "," << genMomPID->at(igen) << "," << genGMomPID->at(igen) << endl;

            if (genMomPID->at(igen) != 23 && genMomPID->at(igen) == 22)
            {
                cout << ev << " " << genPID->at(igen) << " " << genMomPID->at(igen) << " " << genGMomPID->at(igen) << endl;
            }

            // Status Flag requirement: from hard process final state & is prompt final state
            if (((genStatusFlag->at(igen) >> 1) & 1) != 1 || ((genStatusFlag->at(igen) >> 0) & 1) != 1)
                continue;

            if (genMomPID->at(igen) != 23 && genMomPID->at(igen) == 22)
            // if (genMomPID->at(igen) != 23)
            {
                cout << ev << " " << genPID->at(igen) << " " << genMomPID->at(igen) << " " << genGMomPID->at(igen) << endl;
            }

            // if(genMomPID->at(igen) == 23)
            //     vec_Mu.push_back(igen);
            vec_Mu.push_back(igen);

            hM_genMuPt->Fill(genPt->at(igen));
            hM_genMuEta->Fill(genEta->at(igen));
            hM_genMuPhi->Fill(genPhi->at(igen));
        }

        // after imposing Status Flag requirement, do vec_Mu.size() != 2?
        if (vec_Mu.size() != 2)
        {
            for (size_t i = 0; i < vec_Mu.size(); i++)
            {
                cout << "HMMMMM...... vec_Mu size = " << vec_Mu.size() << endl;
                cout << genPID->at(vec_Mu[i]) << " " << genMomPID->at(vec_Mu[i]) << " " << genGMomPID->at(vec_Mu[i]) << endl;
            }
            continue;
        }

        if (genGMomPID->at(vec_Mu[0]) != genGMomPID->at(vec_Mu[1]))
            cout << "genmu1's GMom != genmu2's GMom" << endl;

        fp << ev << "," << genPID->at(vec_Mu[0]) << "," << genMomPID->at(vec_Mu[0]) << "," << genGMomPID->at(vec_Mu[0]) << "," << genPID->at(vec_Mu[1]) << "," << genMomPID->at(vec_Mu[1]) << "," << genGMomPID->at(vec_Mu[1]) << endl;

        // for (size_t i = 0; i < vec_Mu.size(); i++)
        // {
        //     cout << ev << " " << genPID->at(vec_Mu[i]) << " " << genMomPID->at(vec_Mu[i]) << " " << genGMomPID->at(vec_Mu[i]) << endl;

        // }

        TLorentzVector mu1, mu2;
        mu1.SetPtEtaPhiM(genPt->at(vec_Mu[0]), genEta->at(vec_Mu[0]), genPhi->at(vec_Mu[0]), 0.1057);
        mu2.SetPtEtaPhiM(genPt->at(vec_Mu[1]), genEta->at(vec_Mu[1]), genPhi->at(vec_Mu[1]), 0.1057);
        hM_genMmm->Fill((mu1 + mu2).M());
        hM_genMomMass->Fill(genMomMass->at(vec_Mu[0]));
        hM_genMomMass_lowM->Fill(genMomMass->at(vec_Mu[0]));
    }

    // Draw_single1Dhist(TH1F *hist, bool norm1, bool logy, const char *XaxisName, const char *Ytitle_unit, const char *legtext, const char *outname)
    system("mkdir -p plot");
    Draw_single1Dhist(hM_genMuPt, false, true, "pT (GeV)", "GeV", 1, "gen-level pT", "./plot/genPt_PID13_hardproc_isprompt_mZeq2mu");
    Draw_single1Dhist(hM_genMuEta, false, false, "#eta", "", 1, "gen-level Eta", "./plot/genEta_PID13_hardproc_isprompt_mZeq2mu");
    Draw_single1Dhist(hM_genMuPhi, false, false, "#phi", "", 1, "gen-level Phi", "./plot/genPhi_PID13_hardproc_isprompt_mZeq2mu");
    Draw_single1Dhist(hM_genMmm, false, true, "M_{#mu#mu} (GeV)", "GeV", 1, "gen-level di-muon mass", "./plot/genMmm_PID13_MomPID23_hardproc_isprompt_mZeq2mu");
    Draw_single1Dhist(hM_genMomMass, false, true, "M_{#mu#mu} (GeV)", "GeV", 1, "gen-level Mom mass", "./plot/genMomMass_PID13_MomPID23_hardproc_isprompt_mZeq2mu");
    Draw_single1Dhist(hM_genMomMass_lowM, false, true, "M_{#mu#mu} (GeV)", "GeV", 2, "gen-level Mom mass", "./plot/genMomMass_lowMass_PID13_MomPID23_hardproc_isprompt_mZeq2mu");
}