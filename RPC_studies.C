#define RPC_studies_cxx
#include "RPC_studies.h"
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include "TMath.h"


void RPC_studies::Loop()
{

gStyle->SetOptStat(1);
gStyle->SetOptFit(1);

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
	h_residuals_MB1_phi_bending->GetXaxis()->SetTitle("residuals");
	h_residuals_MB2_phi_bending->GetXaxis()->SetTitle("residuals");
	h_MB1_phi_bending_RPC->GetXaxis()->SetTitle("#phiB");
	h_MB2_phi_bending_RPC->GetXaxis()->SetTitle("#phiB");
	h_MB1_phi_bending_DT->GetXaxis()->SetTitle("#phiB");
	h_MB2_phi_bending_DT->GetXaxis()->SetTitle("#phiB");
	h_phi_RPC_vs_DT_MB1->GetXaxis()->SetTitle("DT #phi");
	h_phi_RPC_vs_DT_MB1->GetYaxis()->SetTitle("RPC #phi");
	h_phi_RPC_vs_DT_MB2->GetXaxis()->SetTitle("DT #phi");
	h_phi_RPC_vs_DT_MB2->GetYaxis()->SetTitle("RPC #phi");
	h_phiB_RPC_vs_DT_MB1->GetXaxis()->SetTitle("DT #phiB");
	h_phiB_RPC_vs_DT_MB1->GetYaxis()->SetTitle("RPC #phiB");
	h_phiB_RPC_vs_DT_MB2->GetXaxis()->SetTitle("DT #phiB");
	h_phiB_RPC_vs_DT_MB2->GetYaxis()->SetTitle("RPC #phiB");
	h_MB1_phi_RPC->GetXaxis()->SetTitle("global #phi");
	h_MB2_phi_RPC->GetXaxis()->SetTitle("global #phi");
	h_MB1_phi_DT->GetXaxis()->SetTitle("global #phi");
	h_MB2_phi_DT->GetXaxis()->SetTitle("global #phi");
	h_phi_residuals_MB1->GetXaxis()->SetTitle("residuals");
	h_phi_residuals_MB2->GetXaxis()->SetTitle("residuals");

   if (fChain == 0) return;
   
// 	gStyle->SetOptStat(0);
		
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

				dt_ring = DTextrapolatedOnRPCRing->at(i)+2;
				dt_station = DTextrapolatedOnRPCStation->at(i)-1;
				dt_sector = DTextrapolatedOnRPCSector->at(i)-1;
				
				if(DTextrapolatedOnRPCLayer->at(i) == 1){
					den[0][dt_ring][dt_station][dt_sector]++;
					stripw[0][dt_ring][dt_station][dt_sector] = DTextrapolatedOnRPCStripw->at(i);
					DT_glob_phi[0][dt_ring][dt_station][dt_sector] = DTextrapolatedOnRPCGlobPhi->at(i);
					DT_glob_x[0][dt_ring][dt_station][dt_sector] = DTextrapolatedOnRPCGlobX->at(i);
					DT_glob_y[0][dt_ring][dt_station][dt_sector] = DTextrapolatedOnRPCGlobY->at(i);
				}
		
				if(DTextrapolatedOnRPCLayer->at(i) == 2){
					den[1][dt_ring][dt_station][dt_sector]++;
					stripw[1][dt_ring][dt_station][dt_sector] = DTextrapolatedOnRPCStripw->at(i);
					DT_glob_phi[1][dt_ring][dt_station][dt_sector] = DTextrapolatedOnRPCGlobPhi->at(i);
					DT_glob_x[1][dt_ring][dt_station][dt_sector] = DTextrapolatedOnRPCGlobX->at(i);
					DT_glob_y[1][dt_ring][dt_station][dt_sector] = DTextrapolatedOnRPCGlobY->at(i);

				}
				
				min_layer[0] = 100;
				min_layer[1] = 100;
				
				for(int j = 0; j < NirpcrechitsTwinMux; j++){
					if((RpcRecHitTwinMuxLayer->at(j) == DTextrapolatedOnRPCLayer->at(i)) && 
						(RpcRecHitTwinMuxStation->at(j) == DTextrapolatedOnRPCStation->at(i)) && 
						(RpcRecHitTwinMuxSector->at(j) == DTextrapolatedOnRPCSector->at(i)) && 
						(RpcRecHitTwinMuxRing->at(j) == DTextrapolatedOnRPCRing->at(i))){
				
						if(DTextrapolatedOnRPCLayer->at(i) == 1 && den[0][dt_ring][dt_station][dt_sector] > 0){					
							
							dist_layer[0] = abs(DTextrapolatedOnRPCLocX->at(i) - RpcRechitTwinMuxLocX->at(j));
							
							if(dist_layer[0] < min_layer[0]){
								min_layer[0] = dist_layer[0];
								rpc_ring_layer[0] = RpcRecHitTwinMuxRing->at(j) + 2;
								rpc_station_layer[0] = RpcRecHitTwinMuxStation->at(j) - 1;
								rpc_sector_layer[0] = RpcRecHitTwinMuxSector->at(j) - 1;
								RPC_DT_dist[0][rpc_ring_layer[0]][rpc_station_layer[0]][rpc_sector_layer[0]] = min_layer[0];
								RPC_cluSize[0][rpc_ring_layer[0]][rpc_station_layer[0]][rpc_sector_layer[0]] = RpcRecHitTwinMuxClusterSize->at(j);
								RPC_BX[0][rpc_ring_layer[0]][rpc_station_layer[0]][rpc_sector_layer[0]] = RpcRecHitTwinMuxBx->at(j);
// 								RPC_glob_x_loc_1[rpc_ring_layer[0]][rpc_station_layer[0]][rpc_sector_layer[0]] = DTextrapolatedOnRPCLocX->at(i) - RpcRechitTwinMuxLocX->at(j);
								RPC_glob_phi[0][rpc_ring_layer[0]][rpc_station_layer[0]][rpc_sector_layer[0]] = RpcRechitTwinMuxGlobPhi->at(j);
								RPC_glob_x[0][rpc_ring_layer[0]][rpc_station_layer[0]][rpc_sector_layer[0]] = RpcRechitTwinMuxGlobX->at(j);
								RPC_glob_y[0][rpc_ring_layer[0]][rpc_station_layer[0]][rpc_sector_layer[0]] = RpcRechitTwinMuxGlobY->at(j);
							}
						}
				
						if(DTextrapolatedOnRPCLayer->at(i) == 2 && den[1][dt_ring][dt_station][dt_sector] > 0){			
							
							
							dist_layer[1] = abs(DTextrapolatedOnRPCLocX->at(i) - RpcRechitTwinMuxLocX->at(j));
							if(dist_layer[1] < min_layer[1]){
								min_layer[1] = dist_layer[1];
								rpc_ring_layer[1] = RpcRecHitTwinMuxRing->at(j) + 2;
								rpc_station_layer[1] = RpcRecHitTwinMuxStation->at(j) - 1;
								rpc_sector_layer[1] = RpcRecHitTwinMuxSector->at(j) - 1;
								RPC_DT_dist[1][rpc_ring_layer[1]][rpc_station_layer[1]][rpc_sector_layer[1]] = min_layer[1];
								RPC_cluSize[1][rpc_ring_layer[1]][rpc_station_layer[1]][rpc_sector_layer[1]] = RpcRecHitTwinMuxClusterSize->at(j);
								RPC_BX[1][rpc_ring_layer[1]][rpc_station_layer[1]][rpc_sector_layer[1]] = RpcRecHitTwinMuxBx->at(j);
// 								RPC_glob_x_loc_2[rpc_ring_layer[0]][rpc_station_layer[0]][rpc_sector_layer[0]] = DTextrapolatedOnRPCLocX->at(i) - RpcRechitTwinMuxLocX->at(j);
								RPC_glob_phi[1][rpc_ring_layer[1]][rpc_station_layer[1]][rpc_sector_layer[1]] = RpcRechitTwinMuxGlobPhi->at(j);
								RPC_glob_x[1][rpc_ring_layer[1]][rpc_station_layer[1]][rpc_sector_layer[1]] = RpcRechitTwinMuxGlobX->at(j);
								RPC_glob_y[1][rpc_ring_layer[1]][rpc_station_layer[1]][rpc_sector_layer[1]] = RpcRechitTwinMuxGlobY->at(j);
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
					if(den[0][ii][jj][zz] > 0 && den[1][ii][jj][zz] > 0){  // den: look for one DT_recHit on layer 1 and one on layer 2
						if(jj == 0){
	    	  				h_den_station_1->Fill(zz+1, ii-2);
		    	  			h_dist_MB1_layer_1->Fill(RPC_DT_dist[0][ii][jj][zz]);
  		    				h_dist_MB1_layer_2->Fill(RPC_DT_dist[1][ii][jj][zz]);
  	    					station_1_den++;
  	    					h_den_eta_station_1->Fill(eta_mu);
  	    					h_den_phi_station_1->Fill(phi_mu);
      					}
      					else{				
	      					h_den_station_2->Fill(zz+1, ii-2);
		    	  			h_dist_MB2_layer_1->Fill(RPC_DT_dist[0][ii][jj][zz]);
  		    				h_dist_MB2_layer_2->Fill(RPC_DT_dist[1][ii][jj][zz]);
  	    					station_2_den++;
  	    					h_den_eta_station_2->Fill(eta_mu);
  	    					h_den_phi_station_2->Fill(phi_mu);
  	    				}	      							
      					if(RPC_DT_dist[0][ii][jj][zz] < (rangestrips+RPC_cluSize[0][ii][jj][zz]*0.5)*stripw[0][ii][jj][zz] && RPC_cluSize[0][ii][jj][zz]<=clsCut && // num: looking for one RPC_recHit on layer 1 near DT_recHit &&
     					   RPC_DT_dist[1][ii][jj][zz] < (rangestrips+RPC_cluSize[1][ii][jj][zz]*0.5)*stripw[1][ii][jj][zz] && RPC_cluSize[1][ii][jj][zz]<=clsCut && // num: looking for one RPC_recHit on layer 2 near DT_recHit &&
     					   RPC_BX[0][ii][jj][zz] == RPC_BX[1][ii][jj][zz]){ // RPC_recHit on both layers must have the same BX 

	   						if(jj == 0){
	      						h_num_station_1->Fill(zz+1, ii-2);
			  	  				station_1_num++;
  	    						h_num_eta_station_1->Fill(eta_mu);
  	    						h_num_phi_station_1->Fill(phi_mu);
	  	    					//phi
	  	    					h_MB1_phi_RPC->Fill(RPC_glob_phi[0][ii][jj][zz]);
	  	    					h_MB1_phi_RPC->Fill(RPC_glob_phi[1][ii][jj][zz]);
	  	    					h_MB1_phi_DT->Fill(DT_glob_phi[0][ii][jj][zz]);
	  	    					h_MB1_phi_DT->Fill(DT_glob_phi[1][ii][jj][zz]);
	  	    					h_phi_residuals_MB1->Fill(DT_glob_phi[0][ii][jj][zz] - RPC_glob_phi[0][ii][jj][zz]);
  		    					h_phi_residuals_MB1->Fill(DT_glob_phi[1][ii][jj][zz] - RPC_glob_phi[1][ii][jj][zz]);
  		    					//phi
  		    					//phiB
  		    					DT_phi_bending = TMath::ATan((DT_glob_x[1][ii][jj][zz] - DT_glob_x[0][ii][jj][zz]) / (DT_glob_y[1][ii][jj][zz] - DT_glob_y[0][ii][jj][zz])) - (DT_glob_phi[0][ii][jj][zz] + DT_glob_phi[1][ii][jj][zz])/2;
  		    					RPC_phi_bending = TMath::ATan((RPC_glob_x[1][ii][jj][zz] - RPC_glob_x[0][ii][jj][zz]) / (RPC_glob_y[1][ii][jj][zz] - RPC_glob_y[0][ii][jj][zz])) - (RPC_glob_phi[0][ii][jj][zz]+RPC_glob_phi[1][ii][jj][zz])/2;
  		    					h_MB1_phi_bending_RPC->Fill(RPC_phi_bending);
  		    					h_MB1_phi_bending_DT->Fill(DT_phi_bending);
  		    					h_residuals_MB1_phi_bending->Fill(DT_phi_bending - RPC_phi_bending);
  		    					//phiB
  		    					//TH2
  		    					h_phi_RPC_vs_DT_MB1->Fill(DT_glob_phi[0][ii][jj][zz], RPC_glob_phi[0][ii][jj][zz]);
  		    					h_phi_RPC_vs_DT_MB1->Fill(DT_glob_phi[1][ii][jj][zz], RPC_glob_phi[1][ii][jj][zz]);
  		    					h_phiB_RPC_vs_DT_MB1->Fill(DT_phi_bending, RPC_phi_bending);
  		    					//TH2
	  	    				}
	      					else{
	      						h_num_station_2->Fill(zz+1, ii-2);
	  	    					station_2_num++;
	  	    					h_num_eta_station_2->Fill(eta_mu);
  	    						h_num_phi_station_2->Fill(phi_mu);
	  	    					//phi
	  	    					h_MB2_phi_RPC->Fill(RPC_glob_phi[0][ii][jj][zz]);
	  	    					h_MB2_phi_RPC->Fill(RPC_glob_phi[1][ii][jj][zz]);
	  	    					h_MB2_phi_DT->Fill(DT_glob_phi[0][ii][jj][zz]);
	  	    					h_MB2_phi_DT->Fill(DT_glob_phi[1][ii][jj][zz]);
	  	    					h_phi_residuals_MB2->Fill(DT_glob_phi[0][ii][jj][zz] - RPC_glob_phi[0][ii][jj][zz]);
  		    					h_phi_residuals_MB2->Fill(DT_glob_phi[1][ii][jj][zz] - RPC_glob_phi[1][ii][jj][zz]);
  		    					//phi
  		    					//phiB
  		    					DT_phi_bending = TMath::ATan((DT_glob_x[1][ii][jj][zz] - DT_glob_x[0][ii][jj][zz]) / (DT_glob_y[1][ii][jj][zz] - DT_glob_y[0][ii][jj][zz])) - (DT_glob_phi[0][ii][jj][zz] + DT_glob_phi[1][ii][jj][zz])/2;
  		    					RPC_phi_bending = TMath::ATan((RPC_glob_x[1][ii][jj][zz] - RPC_glob_x[0][ii][jj][zz]) / (RPC_glob_y[1][ii][jj][zz] - RPC_glob_y[0][ii][jj][zz])) - (RPC_glob_phi[0][ii][jj][zz]+RPC_glob_phi[1][ii][jj][zz])/2;
  		    					h_MB2_phi_bending_RPC->Fill(RPC_phi_bending);
  		    					h_MB2_phi_bending_DT->Fill(DT_phi_bending);
  		    					h_residuals_MB2_phi_bending->Fill(DT_phi_bending - RPC_phi_bending);
  		    					//phiB
  		    					//TH2
  		    					h_phi_RPC_vs_DT_MB2->Fill(DT_glob_phi[0][ii][jj][zz], RPC_glob_phi[0][ii][jj][zz]);
  		    					h_phi_RPC_vs_DT_MB2->Fill(DT_glob_phi[1][ii][jj][zz], RPC_glob_phi[1][ii][jj][zz]);
  		    					h_phiB_RPC_vs_DT_MB2->Fill(DT_phi_bending, RPC_phi_bending);
  		    					//TH2
	  	    				}
	      				} // if on RPC
					} //if on DT
					den[0][ii][jj][zz] = 0;
					den[1][ii][jj][zz] = 0;
					RPC_DT_dist[0][ii][jj][zz] = 0;
					RPC_DT_dist[1][ii][jj][zz] = 0;
				} //for on sector
			} //for on ring 
		} // for on station   
   } // loop on event
   
//    std::cout<<DT_segment<<"\t"<<RPC_segment<<"\t"<<RPC_segment_1<<"\t"<<RPC_segment_5<<std::endl;

   std::cout<<station_1_num<<"\t"<<station_1_den<<"\t"<<(float)station_1_num/(float)station_1_den<<std::endl;
   std::cout<<station_2_num<<"\t"<<station_2_den<<"\t"<<(float)station_2_num/(float)station_2_den<<std::endl;

   TCanvas *DT_muon = new TCanvas("DT from muon", "DT from muon", 210,45,750,500);
   h_dist_DT_muon->Draw();

   TCanvas *MB1_layer_1 = new TCanvas("MB1 Layer 1", "MB1 Layer 1", 210,45,750,500);
   h_dist_MB1_layer_1->Draw();

   TCanvas *MB1_layer_2 = new TCanvas("MB1 Layer 2", "MB1 Layer 2", 210,45,750,500);
   h_dist_MB1_layer_2->Draw();
   
   TCanvas *MB2_layer_1 = new TCanvas("MB2 Layer 1", "MB2 Layer 1", 210,45,750,500);
   h_dist_MB2_layer_1->Draw();

   TCanvas *MB2_layer_2 = new TCanvas("MB2 Layer 2", "MB2 Layer 2", 210,45,750,500);
   h_dist_MB2_layer_2->Draw();
   
   TCanvas *num_station_1 = new TCanvas("num station 1", "num station 1", 210,45,750,500);
   h_num_station_1->Draw("COLZ"); 
   h_num_station_1->Draw("TEXT SAME");   
   h_num_station_1->SetMinimum(0);
   h_num_station_1->SetMaximum(50);  

   TCanvas *num_station_2 = new TCanvas("num station 2", "num station 2", 210,45,750,500);
   h_num_station_2->Draw("COLZ");
   h_num_station_2->Draw("TEXT SAME");   
   h_num_station_2->SetMinimum(0);
   h_num_station_2->SetMaximum(50);
   
   TCanvas *den_station_1 = new TCanvas("den station 1", "den station 1", 210,45,750,500);
   h_den_station_1->Draw("COLZ"); 
   h_den_station_1->Draw("TEXT SAME");     
   h_den_station_1->SetMinimum(0);
   h_den_station_1->SetMaximum(50);  

   TCanvas *den_station_2 = new TCanvas("den station 2", "den station 2", 210,45,750,500);
   h_den_station_2->Draw("COLZ");
   h_den_station_2->Draw("TEXT SAME");     
   h_den_station_2->SetMinimum(0);
   h_den_station_2->SetMaximum(50);
  
   if(TEfficiency::CheckConsistency(*h_num_station_1, *h_den_station_1))   
			h_eff_station_1 = new TEfficiency(*h_num_station_1, *h_den_station_1);
   if(TEfficiency::CheckConsistency(*h_num_station_2, *h_den_station_2))   
		    h_eff_station_2 = new TEfficiency(*h_num_station_2, *h_den_station_2);
   
   TCanvas *eff_station_1 = new TCanvas("eff station 1", "eff station 1", 210,45,750,500);
   h_eff_station_1->Draw("COLZ"); 
   h_eff_station_1->Draw("TEXT SAME");  
   eff_station_1->Print("Eff_MB1.png");
   eff_station_1->Print("Eff_MB1.pdf");

   TCanvas *eff_station_2 = new TCanvas("eff station 2", "eff station 2", 210,45,750,500);
   h_eff_station_2->Draw("COLZ");
   h_eff_station_2->Draw("TEXT SAME");  
   eff_station_2->Print("Eff_MB2.png");
   eff_station_2->Print("Eff_MB2.pdf");


   if(TEfficiency::CheckConsistency(*h_num_eta_station_1, *h_den_eta_station_1))   
			h_eff_eta_station_1 = new TEfficiency(*h_num_eta_station_1, *h_den_eta_station_1);
   if(TEfficiency::CheckConsistency(*h_num_eta_station_2, *h_den_eta_station_2))   
		    h_eff_eta_station_2 = new TEfficiency(*h_num_eta_station_2, *h_den_eta_station_2);
   
   TCanvas *eff_eta = new TCanvas("eff vs #eta", "eff vs #eta", 210,45,750,500);
   h_eff_eta_station_1->SetLineColor(kRed);
   h_eff_eta_station_2->SetLineColor(kBlue);
   h_eff_eta_station_1->Draw(); 
   h_eff_eta_station_2->Draw("same");   
   TLegend *legend_eff_eta = new TLegend(0.85,0.8,0.9,0.9);
   legend_eff_eta->AddEntry(h_eff_eta_station_1, "MB1", "lep");
   legend_eff_eta->AddEntry(h_eff_eta_station_2,"MB2", "lep");	
   legend_eff_eta->Draw();
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
   h_eff_phi_station_2->Draw("same");
   TLegend *legend_eff_phi = new TLegend(0.85,0.8,0.9,0.9);
   legend_eff_phi->AddEntry(h_eff_phi_station_1, "MB1", "lep");
   legend_eff_phi->AddEntry(h_eff_phi_station_2,"MB2", "lep");	
   legend_eff_phi->Draw();
   eff_phi->Print("Eff_vs_phi.png");
   eff_phi->Print("Eff_vs_phi.pdf");


   TCanvas *MB1_phi_RPC = new TCanvas("#phi on MB1 RPC", "#phi on MB1 RPC", 210,45,750,500);
   h_MB1_phi_RPC->Draw();
   MB1_phi_RPC->Print("phi_RPC_MB1.pdf");
   
   TCanvas *MB2_phi_RPC = new TCanvas("#phi on MB2 RPC", "#phi on MB2 RPC", 210,45,750,500);
   h_MB2_phi_RPC->Draw();
   MB2_phi_RPC->Print("phi_RPC_MB2.pdf");
   
   TCanvas *MB1_phi_DT = new TCanvas("#phi on MB1 DT", "#phi on MB1 DT", 210,45,750,500);
   h_MB1_phi_DT->Draw();
   MB1_phi_DT->Print("phi_DT_MB1.pdf");
   
   TCanvas *MB2_phi_DT = new TCanvas("#phi on MB2 DT", "#phi on MB2 DT", 210,45,750,500);
   h_MB2_phi_DT->Draw();
   MB2_phi_DT->Print("phi_DT_MB2.pdf");   
   
   TCanvas *phi_residuals_MB1 = new TCanvas("MB1 layer", "MB1 layer", 210,45,750,500);
   fit_min = h_phi_residuals_MB1->GetMean() - h_phi_residuals_MB1->GetRMS();
   fit_max = h_phi_residuals_MB1->GetMean() + h_phi_residuals_MB1->GetRMS();
   TF1 *f1 = new TF1("f1","gaus",fit_min, fit_max);
   h_phi_residuals_MB1->Fit("f1", "R");
   h_phi_residuals_MB1->Draw();
   phi_residuals_MB1->Print("phi_residuals_MB1.pdf");
   
   TCanvas *phi_residuals_MB2 = new TCanvas("MB2 layer", "MB2 layer", 210,45,750,500);
   fit_min = h_phi_residuals_MB2->GetMean() - h_phi_residuals_MB2->GetRMS();
   fit_max = h_phi_residuals_MB2->GetMean() + h_phi_residuals_MB2->GetRMS();
   TF1 *f2 = new TF1("f2","gaus",fit_min, fit_max);
   h_phi_residuals_MB2->Fit("f2", "R");
   h_phi_residuals_MB2->Draw();
   phi_residuals_MB2->Print("phi_residuals_MB2.pdf");   

   TCanvas *phi_bending_MB1_DT = new TCanvas("#phi bending MB1 DT", "#phi bending MB1 DT", 210,45,750,500);
   h_MB1_phi_bending_DT->Draw();
   phi_bending_MB1_DT->Print("phi_bending_MB1_DT.pdf");

   TCanvas *phi_bending_MB2_DT = new TCanvas("#phi bending MB2 DT", "#phi bending MB2 DT", 210,45,750,500);
   h_MB2_phi_bending_DT->Draw();
   phi_bending_MB2_DT->Print("phi_bending_MB2_DT.pdf");
   
   TCanvas *phi_bending_MB1_RPC = new TCanvas("#phi bending MB1 RPC", "#phi bending MB1 RPC", 210,45,750,500);
   h_MB1_phi_bending_RPC->Draw();
   phi_bending_MB1_RPC->Print("phi_bending_MB1_RPC.pdf");

   TCanvas *phi_bending_MB2_RPC = new TCanvas("#phi bending MB2 RPC", "#phi bending MB2 RPC", 210,45,750,500);
   h_MB2_phi_bending_RPC->Draw();
   phi_bending_MB2_RPC->Print("phi_bending_MB2_RPC.pdf");
   
   TCanvas *phi_bending_residuals_MB1 = new TCanvas("#phi bending residuals MB1", "#phi bending residuals MB1", 210,45,750,500);
   fit_min = h_residuals_MB1_phi_bending->GetMean() - h_residuals_MB1_phi_bending->GetRMS();
   fit_max = h_residuals_MB1_phi_bending->GetMean() + h_residuals_MB1_phi_bending->GetRMS();
   TF1 *f3 = new TF1("f3","gaus",fit_min, fit_max);
   h_residuals_MB1_phi_bending->Fit("f3", "R");
   h_residuals_MB1_phi_bending->Draw();
   phi_bending_residuals_MB1->Print("phi_bending_residuals_MB1.pdf");

   TCanvas *phi_bending_residuals_MB2 = new TCanvas("#phi bending residuals MB2", "#phi bending residuals MB2", 210,45,750,500);
   fit_min = h_residuals_MB2_phi_bending->GetMean() - h_residuals_MB2_phi_bending->GetRMS();
   fit_max = h_residuals_MB2_phi_bending->GetMean() + h_residuals_MB2_phi_bending->GetRMS();
   TF1 *f4 = new TF1("f4","gaus",fit_min, fit_max);
   h_residuals_MB2_phi_bending->Fit("f4", "R");
   h_residuals_MB2_phi_bending->Draw();
   phi_bending_residuals_MB2->Print("phi_bending_residuals_MB2.pdf");

   
   TCanvas *RPC_vs_DT_phi_MB1 = new TCanvas("#phi: RPC vs DT MB1", "#phi: RPC vs DT MB1", 210,45,750,500);
   h_phi_RPC_vs_DT_MB1->Draw();
   RPC_vs_DT_phi_MB1->Print("RPC_vs_DT_phi_MB1.pdf");
   
   TCanvas *RPC_vs_DT_phi_MB2 = new TCanvas("#phi: RPC vs DT MB2", "#phi: RPC vs DT MB2", 210,45,750,500);
   h_phi_RPC_vs_DT_MB2->Draw();
   RPC_vs_DT_phi_MB2->Print("RPC_vs_DT_phi_MB2.pdf");
   
   TCanvas *RPC_vs_DT_phiB_MB1 = new TCanvas("#phiB: RPC vs DT MB1", "#phiB: RPC vs DT MB1", 210,45,750,500);
   h_phiB_RPC_vs_DT_MB1->Draw();
   RPC_vs_DT_phiB_MB1->Print("RPC_vs_DT_phiB_MB1.pdf");
   
   TCanvas *RPC_vs_DT_phiB_MB2 = new TCanvas("#phiB: RPC vs DT MB2", "#phiB: RPC vs DT MB2", 210,45,750,500);
   h_phiB_RPC_vs_DT_MB2->Draw();
   RPC_vs_DT_phiB_MB2->Print("RPC_vs_DT_phiB_MB2.pdf");



std::cout<<count<<std::endl;
}
