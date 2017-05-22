#define RPC_studies_cxx
#include "RPC_studies.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>


void RPC_studies::Loop()
{
//   In a ROOT session, you can do:
//      root> .L RPC_studies.C
//      root> RPC_studies t
//      root> t.GetEntry(12); // Fill t data members with entry number 12
//      root> t.Show();       // Show values of entry 12
//      root> t.Show(16);     // Read and show values of entry 16
//      root> t.Loop();       // Loop on all entries
//

//     This is the loop skeleton where:
//    jentry is the global entry number in the chain
//    ientry is the entry number in the current Tree
//  Note that the argument to GetEntry must be:
//    jentry for TChain::GetEntry
//    ientry for TTree::GetEntry and TBranch::GetEntry
//
//       To read only selected branches, Insert statements like:
// METHOD1:
//    fChain->SetBranchStatus("*",0);  // disable all branches
//    fChain->SetBranchStatus("branchname",1);  // activate branchname
// METHOD2: replace line
//    fChain->GetEntry(jentry);       //read all branches
//by  b_branchname->GetEntry(ientry); //read only this branch
	h_dist_DT_muon->GetXaxis()->SetTitle("#DeltaR");
	h_dist_MB1_layer_1->GetXaxis()->SetTitle("#DeltaR");
	h_dist_MB1_layer_2->GetXaxis()->SetTitle("#DeltaR");
	h_dist_MB2_layer_1->GetXaxis()->SetTitle("#DeltaR");
	h_dist_MB2_layer_2->GetXaxis()->SetTitle("#DeltaR");
	h_num_station_1->GetYaxis()->SetTitle("Ring/Wheel");
	h_num_station_2->GetYaxis()->SetTitle("Ring/Wheel");
	h_den_station_1->GetYaxis()->SetTitle("Ring/Wheel");
	h_den_station_2->GetYaxis()->SetTitle("Ring/Wheel");
	h_num_station_1->GetXaxis()->SetTitle("Sector");
	h_num_station_2->GetXaxis()->SetTitle("Sector");
	h_den_station_1->GetXaxis()->SetTitle("Sector");
	h_den_station_2->GetXaxis()->SetTitle("Sector");
	h_num_eta_station_1->GetXaxis()->SetTitle("#eta");
	h_num_eta_station_2->GetXaxis()->SetTitle("#eta");
	h_den_eta_station_1->GetXaxis()->SetTitle("#eta");
	h_den_eta_station_2->GetXaxis()->SetTitle("#eta");
	h_num_eta_station_1->GetYaxis()->SetTitle("efficiency");
	h_num_eta_station_2->GetYaxis()->SetTitle("efficiency");
	h_den_eta_station_1->GetYaxis()->SetTitle("efficiency");
	h_den_eta_station_2->GetYaxis()->SetTitle("efficiency");	
	h_num_phi_station_1->GetXaxis()->SetTitle("#phi");
	h_num_phi_station_2->GetXaxis()->SetTitle("#phi");
	h_den_phi_station_1->GetXaxis()->SetTitle("#phi");
	h_den_phi_station_2->GetXaxis()->SetTitle("#phi");
	h_num_phi_station_1->GetYaxis()->SetTitle("efficiency");
	h_num_phi_station_2->GetYaxis()->SetTitle("efficiency");
	h_den_phi_station_1->GetYaxis()->SetTitle("efficiency");
	h_den_phi_station_2->GetYaxis()->SetTitle("efficiency");

	

   if (fChain == 0) return;
   
	gStyle->SetOptStat(0);
		
   Long64_t nentries = fChain->GetEntriesFast();

   Long64_t nbytes = 0, nb = 0;
   
   std::cout<<nentries<<std::endl;
   
   for (Long64_t jentry=0; jentry<nentries;jentry++) {
      Long64_t ientry = LoadTree(jentry);
      if (ientry < 0) break;
      nb = fChain->GetEntry(jentry);   nbytes += nb;
      // if (Cut(ientry) < 0) continue;
      
      if(jentry % 10000 == 0) std::cout<<jentry<<std::endl;
      
      for(int y = 0; y < Ndtsegments; y++){
      		
      		min_DT_muon = 20;
      		match_DT_muon = false;
      		      		
      		for(int h = 0; h < Nmuons; h++){
      		
      			Mu_pt = sqrt(Mu_px->at(h)*Mu_px->at(h) + Mu_py->at(h)*Mu_py->at(h));
      			if(Mu_pt < 5) continue;
				
				delta_phi = abs(dtsegm4D_phi->at(y) - Mu_phi->at(h));
				if(delta_phi > pig) delta_phi=2*pig-delta_phi;
				
				delta_eta = abs(dtsegm4D_eta->at(y)-Mu_eta->at(h));
				
				deltaR = sqrt(delta_phi*delta_phi + delta_eta*delta_eta);
				if(deltaR < min_DT_muon){
					min_DT_muon = deltaR;
					if(Mu_isMuGlobal->at(h)==1 && Mu_isMuTracker->at(h)==1 &&
						fabs(Mu_dxy_glb->at(h))<D0_cut && fabs(Mu_dz_glb->at(h))<Dz_cut &&
						fabs(Mu_eta->at(h))<eta_cut && (Mu_normchi2_glb->at(h)) < muchi2_cut &&
						Mu_pt>pt_cut && Mu_pt<pt_max && 
						Mu_numberOfHits_sta->at(h)>=N_hits_cut && 
						Mu_numberOfPixelHits_glb->at(h)>=npix_cut &&
						Mu_numberOfTrackerHits_glb->at(h)>=ntkr_cut ) {
																		eta_mu = Mu_eta->at(h);
																		phi_mu = Mu_phi->at(h);

					}
				}
      		} // for on RecoMuons --- selected as StandAlone in the ntuples
      		
      		h_dist_DT_muon->Fill(min_DT_muon);
      		
      		if(min_DT_muon < 0.4)
      			match_DT_muon = true;
      			
      		if(!match_DT_muon) continue; // if the DT segment is to far from the muon, it is excluded
      		
      		for(int i = 0; i < NDTsegmentonRPC; i++){
      		
				if((dtsegm4D_station->at(y) != DTextrapolatedOnRPCStation->at(i)) || 
						(dtsegm4D_sector->at(y) != DTextrapolatedOnRPCSector->at(i)) || 
						(dtsegm4D_wheel->at(y) != DTextrapolatedOnRPCRing->at(i))) continue;
      		
      			count++;
      			
      			if(DTextrapolatedOnRPCRegion->at(i) != 0) std::cout<<" DT EXTRAPOLATION OUT OF BARREL !!! PLEASE CHECK "<<std::endl;
      			
      			if(DTextrapolatedOnRPCStation->at(i) == 3){
//       				std::cout<<"MB3: we don't care about it"<<std::endl;
      				continue;
      			}
      			
				if(DTextrapolatedOnRPCLayer->at(i) == 1){
					dt_ring_1 = DTextrapolatedOnRPCRing->at(i)+2;
					dt_station_1 = DTextrapolatedOnRPCStation->at(i)-1;
					dt_sector_1 = DTextrapolatedOnRPCSector->at(i)-1;

					den_1[dt_ring_1][dt_station_1][dt_sector_1]++;
					stripw_1[rpc_ring_1][rpc_station_1][rpc_sector_1] = DTextrapolatedOnRPCStripw->at(i);
				}
		
				if(DTextrapolatedOnRPCLayer->at(i) == 2){
					dt_ring_2 = DTextrapolatedOnRPCRing->at(i)+2;
					dt_station_2 = DTextrapolatedOnRPCStation->at(i)-1;
					dt_sector_2 = DTextrapolatedOnRPCSector->at(i)-1;
				
					den_2[dt_ring_2][dt_station_2][dt_sector_2]++;
					stripw_2[rpc_ring_2][rpc_station_2][rpc_sector_2] = DTextrapolatedOnRPCStripw->at(i);
				}
				
				min_layer_1 = 100;
				min_layer_2 = 100;
				
				for(int j = 0; j < NirpcrechitsTwinMux; j++){
					if((RpcRecHitTwinMuxLayer->at(j) == DTextrapolatedOnRPCLayer->at(i)) && 
						(RpcRecHitTwinMuxStation->at(j) == DTextrapolatedOnRPCStation->at(i)) && 
						(RpcRecHitTwinMuxSector->at(j) == DTextrapolatedOnRPCSector->at(i)) && 
						(RpcRecHitTwinMuxRing->at(j) == DTextrapolatedOnRPCRing->at(i))){
				
						if(DTextrapolatedOnRPCLayer->at(i) == 1 && den_1[dt_ring_1][dt_station_1][dt_sector_1] > 0){					
							dist_layer_1 = abs(DTextrapolatedOnRPCLocX->at(i) - RpcRechitTwinMuxLocX->at(j));
							if(dist_layer_1 < min_layer_1){
								min_layer_1 = dist_layer_1;
								rpc_ring_1 = RpcRecHitTwinMuxRing->at(j) + 2;
								rpc_station_1 = RpcRecHitTwinMuxStation->at(j) - 1;
								rpc_sector_1 = RpcRecHitTwinMuxSector->at(j) - 1;
								dist_1[rpc_ring_1][rpc_station_1][rpc_sector_1] = min_layer_1;
								cluSize_1[rpc_ring_1][rpc_station_1][rpc_sector_1] = RpcRecHitTwinMuxClusterSize->at(j);
								BX_1[rpc_ring_1][rpc_station_1][rpc_sector_1] = RpcRecHitTwinMuxBx->at(j);
		
							}
						}
				
						if(DTextrapolatedOnRPCLayer->at(i) == 2 && den_2[dt_ring_2][dt_station_2][dt_sector_2] > 0){			
							dist_layer_2 = abs(DTextrapolatedOnRPCLocX->at(i) - RpcRechitTwinMuxLocX->at(j));
							if(dist_layer_2 < min_layer_2){
								min_layer_2 = dist_layer_2;
								rpc_ring_2 = RpcRecHitTwinMuxRing->at(j) + 2;
								rpc_station_2 = RpcRecHitTwinMuxStation->at(j) - 1;
								rpc_sector_2 = RpcRecHitTwinMuxSector->at(j) - 1;
								dist_2[rpc_ring_2][rpc_station_2][rpc_sector_2] = min_layer_2;
								cluSize_2[rpc_ring_2][rpc_station_2][rpc_sector_2] = RpcRecHitTwinMuxClusterSize->at(j);
								BX_2[rpc_ring_2][rpc_station_2][rpc_sector_2] = RpcRecHitTwinMuxBx->at(j);
							}
						}	
				
					} // if DT = RPC
				} // for on RPC rechit TwinMux  
				
// 				std::cout<<" -------------------------------------------------------------------------------------------- "<<std::endl;
			} // for on DT recHit extrapolated in RPC layers 
		} //for on DT segment
		
		for(int jj = 0; jj < 2; jj++){   
			for(int ii = 0; ii < 5; ii++){
				for(int zz = 0; zz < 12; zz++){
					if(den_1[ii][jj][zz] > 0 && den_2[ii][jj][zz] > 0){  // den: look for one DT_recHit on layer 1 and one on layer 2
						if(jj == 0){
	    	  				h_den_station_1->Fill(zz+1, ii-2);
		    	  			h_dist_MB1_layer_1->Fill(dist_1[ii][jj][zz]);
  		    				h_dist_MB1_layer_2->Fill(dist_2[ii][jj][zz]);
  	    					station_1_den++;
  	    					h_den_eta_station_1->Fill(eta_mu);
  	    					h_den_phi_station_1->Fill(phi_mu);
      					}
      					else{				
	      					h_den_station_2->Fill(zz+1, ii-2);
		    	  			h_dist_MB2_layer_1->Fill(dist_1[ii][jj][zz]);
  		    				h_dist_MB2_layer_2->Fill(dist_2[ii][jj][zz]);
  	    					station_2_den++;
  	    					h_den_eta_station_2->Fill(eta_mu);
  	    					h_den_phi_station_2->Fill(phi_mu);
  	    				}	      							
      					if(dist_1[ii][jj][zz] < (rangestrips+cluSize_1[ii][jj][zz]*0.5)*stripw_1[ii][jj][zz] && cluSize_1[ii][jj][zz]<=clsCut && 
      					   dist_2[ii][jj][zz] < (rangestrips+cluSize_2[ii][jj][zz]*0.5)*stripw_2[ii][jj][zz] && cluSize_2[ii][jj][zz]<=clsCut &&
      					   BX_1[ii][jj][zz] == BX_2[ii][jj][zz]){  // num: look for one RPC_recHit on layer 1 near DT_recHit and one on layer 2 near DT_recHit
//       					if(dist_1[ii][jj][zz] < deltaR && dist_2[ii][jj][zz] < deltaR){ // num: look for one RPC_recHit on layer 1 near DT_recHit and one on layer 2 near DT_recHit
	   						if(jj == 0){
	      						h_num_station_1->Fill(zz+1, ii-2);
			  	  				station_1_num++;
  	    						h_num_eta_station_1->Fill(eta_mu);
  	    						h_num_phi_station_1->Fill(phi_mu);
	  	    				}
	      					else{
	      						h_num_station_2->Fill(zz+1, ii-2);
	  	    					station_2_num++;
	  	    					h_num_eta_station_2->Fill(eta_mu);
  	    						h_num_phi_station_2->Fill(phi_mu);
	  	    				}
	      				} // if on RPC
					} //if on DT
					den_1[ii][jj][zz] = 0;
					den_2[ii][jj][zz] = 0;
					dist_1[ii][jj][zz] = 0;
					dist_2[ii][jj][zz] = 0;
				} //for on sector
			} //for on ring 
		} // for on station   
   } // loop on event
   
//    std::cout<<DT_segment<<"\t"<<RPC_segment<<"\t"<<RPC_segment_1<<"\t"<<RPC_segment_5<<std::endl;

   std::cout<<station_1_num<<"\t"<<station_1_den<<"\t"<<(float)station_1_num/(float)station_1_den<<std::endl;
   std::cout<<station_2_num<<"\t"<<station_2_den<<"\t"<<(float)station_2_num/(float)station_2_den<<std::endl;
   
//    TCanvas *DT_muon = new TCanvas("DT from muon", "DT from muon", 210,45,750,500);
//    h_dist_DT_muon->Draw();
// 
//    TCanvas *MB1_layer_1 = new TCanvas("MB1 Layer 1", "MB1 Layer 1", 210,45,750,500);
//    h_dist_MB1_layer_1->Draw();
// 
//    TCanvas *MB1_layer_2 = new TCanvas("MB1 Layer 2", "MB1 Layer 2", 210,45,750,500);
//    h_dist_MB1_layer_2->Draw();
//    
//    TCanvas *MB2_layer_1 = new TCanvas("MB2 Layer 1", "MB2 Layer 1", 210,45,750,500);
//    h_dist_MB2_layer_1->Draw();
// 
//    TCanvas *MB2_layer_2 = new TCanvas("MB2 Layer 2", "MB2 Layer 2", 210,45,750,500);
//    h_dist_MB2_layer_2->Draw();
//    
//    TCanvas *num_station_1 = new TCanvas("num station 1", "num station 1", 210,45,750,500);
//    h_num_station_1->Draw("COLZ"); 
//    h_num_station_1->Draw("TEXT SAME");   
//    h_num_station_1->SetMinimum(0);
//    h_num_station_1->SetMaximum(50);  
// 
//    TCanvas *num_station_2 = new TCanvas("num station 2", "num station 2", 210,45,750,500);
//    h_num_station_2->Draw("COLZ");
//    h_num_station_2->Draw("TEXT SAME");   
//    h_num_station_2->SetMinimum(0);
//    h_num_station_2->SetMaximum(50);
//    
//    TCanvas *den_station_1 = new TCanvas("den station 1", "den station 1", 210,45,750,500);
//    h_den_station_1->Draw("COLZ"); 
//    h_den_station_1->Draw("TEXT SAME");     
//    h_den_station_1->SetMinimum(0);
//    h_den_station_1->SetMaximum(50);  
// 
//    TCanvas *den_station_2 = new TCanvas("den station 2", "den station 2", 210,45,750,500);
//    h_den_station_2->Draw("COLZ");
//    h_den_station_2->Draw("TEXT SAME");     
//    h_den_station_2->SetMinimum(0);
//    h_den_station_2->SetMaximum(50);
//    
//    if(TEfficiency::CheckConsistency(*h_num_station_1, *h_den_station_1))   
// 			h_eff_station_1 = new TEfficiency(*h_num_station_1, *h_den_station_1);
//    if(TEfficiency::CheckConsistency(*h_num_station_2, *h_den_station_2))   
// 		    h_eff_station_2 = new TEfficiency(*h_num_station_2, *h_den_station_2);
//    
//    TCanvas *eff_station_1 = new TCanvas("eff station 1", "eff station 1", 210,45,750,500);
//    h_eff_station_1->Draw("COLZ"); 
//    h_eff_station_1->Draw("TEXT SAME");  
//    eff_station_1->Print("Eff_MB1.png");
//    eff_station_1->Print("Eff_MB1.pdf");
// 
//    TCanvas *eff_station_2 = new TCanvas("eff station 2", "eff station 2", 210,45,750,500);
//    h_eff_station_2->Draw("COLZ");
//    h_eff_station_2->Draw("TEXT SAME");  
//    eff_station_2->Print("Eff_MB2.png");
//    eff_station_2->Print("Eff_MB2.pdf");


   if(TEfficiency::CheckConsistency(*h_num_eta_station_1, *h_den_eta_station_1))   
			h_eff_eta_station_1 = new TEfficiency(*h_num_eta_station_1, *h_den_eta_station_1);
   if(TEfficiency::CheckConsistency(*h_num_eta_station_2, *h_den_eta_station_2))   
		    h_eff_eta_station_2 = new TEfficiency(*h_num_eta_station_2, *h_den_eta_station_2);
   
   TCanvas *eff_eta = new TCanvas("eff vs #eta", "eff vs #eta", 210,45,750,500);
   h_eff_eta_station_1->SetLineColor(kRed);
   h_eff_eta_station_2->SetLineColor(kBlue);
   h_eff_eta_station_1->Draw(); 
//    h_eff_eta_station_2->Draw("same");   
   TLegend *legend_eff_eta = new TLegend(0.85,0.8,0.9,0.9);
   legend_eff_eta->AddEntry(h_eff_eta_station_1, "MB1", "lep");
   legend_eff_eta->AddEntry(h_eff_eta_station_2,"MB2", "lep");	
//    legend_eff_eta->Draw();
   eff_eta->Print("Eff_vs_eta.png");
   eff_eta->Print("Eff_vs_eta.pdf");
   
   if(TEfficiency::CheckConsistency(*h_num_phi_station_1, *h_den_phi_station_1))   
			h_eff_phi_station_1 = new TEfficiency(*h_num_phi_station_1, *h_den_phi_station_1);
   if(TEfficiency::CheckConsistency(*h_num_phi_station_2, *h_den_phi_station_2))   
		    h_eff_phi_station_2 = new TEfficiency(*h_num_phi_station_2, *h_den_phi_station_2);
   
   TCanvas *eff_phi = new TCanvas("eff vs #phi", "eff vs #phi", 210,45,750,500);
   h_eff_phi_station_1->SetLineColor(kRed);
   h_eff_phi_station_2->SetLineColor(kBlue);
   h_eff_phi_station_1->Draw(); 
//    h_eff_phi_station_2->Draw("same");
   TLegend *legend_eff_phi = new TLegend(0.85,0.8,0.9,0.9);
   legend_eff_phi->AddEntry(h_eff_phi_station_1, "MB1", "lep");
   legend_eff_phi->AddEntry(h_eff_phi_station_2,"MB2", "lep");	
//    legend_eff_phi->Draw();
   eff_phi->Print("Eff_vs_phi.png");
   eff_phi->Print("Eff_vs_phi.pdf");



std::cout<<count<<std::endl;
}
