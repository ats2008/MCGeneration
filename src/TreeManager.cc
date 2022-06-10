void TreeManager::ClearAllVar()
{
        for (std::pair<TString, Float_t> element : storageFloat) {
        storageFloat[element.first]=-1e9;
    }

}
void TreeManager::AddParticleBranch(TTree * aTree, TString tag)
{
    
    storageFloat[tag+"_mass"	    ] = 0.0;
    storageFloat[tag+"_isFinal"	    ] = 0.0;
    storageFloat[tag+"_charge"	    ] = 0.0;
    storageFloat[tag+"_pdgid"	    ] = 0.0;
    storageFloat[tag+"_px"	        ] = 0.0;
    storageFloat[tag+"_py"	        ] = 0.0;
    storageFloat[tag+"_pz"	        ] = 0.0;
    storageFloat[tag+"_e"	        ] = 0.0;
    storageFloat[tag+"_pt"	        ] = 0.0;
    storageFloat[tag+"_rapi"        ] = 0.0;
    storageFloat[tag+"_eta"	        ] = 0.0;
    storageFloat[tag+"_phi"	        ] = 0.0;
    storageFloat[tag+"_vProdX"	    ] = 0.0;
    storageFloat[tag+"_vProdY"	    ] = 0.0;
    storageFloat[tag+"_vProdZ"	    ] = 0.0;
    storageFloat[tag+"_vDecayX"	] = 0.0;
    storageFloat[tag+"_vDecayY"	] = 0.0;
    storageFloat[tag+"_vDecayZ"	] = 0.0;
    storageFloat[tag+"_vProdLXY"    ] = 0.0;
    storageFloat[tag+"_vProdR"      ] = 0.0;
    storageFloat[tag+"_vDecayLXY"  ] = 0.0;
    storageFloat[tag+"_vDecayR"    ] = 0.0;
    storageFloat[tag+"_decayDLXY"  ] = 0.0;
    storageFloat[tag+"_decayDR"    ] = 0.0;
    storageFloat[tag+"_decayDZ"    ] = 0.0;
    
    aTree->Branch(tag+"_mass"	   ,  &   storageFloat[tag+"_mass"	    ]  ) ; 
    aTree->Branch(tag+"_isFinal"   ,  &   storageFloat[tag+"_isFinal"	    ]  ) ;         
    aTree->Branch(tag+"_charge"	   ,  &   storageFloat[tag+"_charge"	    ]  ) ;     
    aTree->Branch(tag+"_pdgid"	   ,  &   storageFloat[tag+"_pdgid"	    ]  ) ;     
    aTree->Branch(tag+"_px"	       ,  &   storageFloat[tag+"_px"	        ]  ) ;     
    aTree->Branch(tag+"_py"	       ,  &   storageFloat[tag+"_py"	        ]  ) ;     
    aTree->Branch(tag+"_pz"	       ,  &   storageFloat[tag+"_pz"	        ]  ) ;     
    aTree->Branch(tag+"_e"	       ,  &   storageFloat[tag+"_e"	        ]  ) ;     
    aTree->Branch(tag+"_pt"	       ,  &   storageFloat[tag+"_pt"	        ]  ) ;     
    aTree->Branch(tag+"_rapi"	   ,  &   storageFloat[tag+"_rapi"	        ]  ) ;         
    aTree->Branch(tag+"_eta"	   ,  &   storageFloat[tag+"_eta"	        ]  ) ;         
    aTree->Branch(tag+"_phi"	   ,  &   storageFloat[tag+"_phi"	        ]  ) ;         
    aTree->Branch(tag+"_vProdX"	   ,  &   storageFloat[tag+"_vProdX"	    ]  ) ;     
    aTree->Branch(tag+"_vProdY"	   ,  &   storageFloat[tag+"_vProdY"	    ]  ) ;     
    aTree->Branch(tag+"_vProdZ"	   ,  &   storageFloat[tag+"_vProdZ"	    ]  ) ;     
    aTree->Branch(tag+"_vDecayX"  ,  &   storageFloat[tag+"_vDecayX"	]  ) ; 
    aTree->Branch(tag+"_vDecayY"  ,  &   storageFloat[tag+"_vDecayY"	]  ) ; 
    aTree->Branch(tag+"_vDecayZ"  ,  &   storageFloat[tag+"_vDecayZ"	]  ) ; 
    aTree->Branch(tag+"_vProdLXY"  ,  &   storageFloat[tag+"_vProdLXY"    ]  ) ;     
    aTree->Branch(tag+"_vProdR"    ,  &   storageFloat[tag+"_vProdR"      ]  ) ;     
    aTree->Branch(tag+"_vDecayLXY",  &   storageFloat[tag+"_vDecayLXY"  ]  ) ;    
    aTree->Branch(tag+"_vDecayR"  ,  &   storageFloat[tag+"_vDecayR"    ]  ) ;     
    aTree->Branch(tag+"_decayDLXY",  &   storageFloat[tag+"_decayDLXY"]  ) ;    
    aTree->Branch(tag+"_decayDR"   ,  &   storageFloat[tag+"_decayDR"  ]  ) ;     
    aTree->Branch(tag+"_decayDZ"   ,  &   storageFloat[tag+"_decayDZ" ]  ) ;     

}

void TreeManager::FillParticleBranches(TString tag,Event &event,Int_t idx)
{

    storageFloat[tag+"_mass"	    ] = event[idx].m();
    storageFloat[tag+"_isFinal"	    ] = event[idx].isFinal();
    storageFloat[tag+"_charge"	    ] = event[idx].charge();
    storageFloat[tag+"_pdgid"	    ] = event[idx].id();
    storageFloat[tag+"_px"	        ] = event[idx].px();
    storageFloat[tag+"_py"	        ] = event[idx].py();
    storageFloat[tag+"_pz"	        ] = event[idx].pz();
    storageFloat[tag+"_e"	        ] = event[idx].e();
    storageFloat[tag+"_pt"	        ] = event[idx].pT();
    storageFloat[tag+"_rapi"        ] = event[idx].y();
    storageFloat[tag+"_eta"	        ] = event[idx].eta();
    storageFloat[tag+"_phi"	        ] = event[idx].phi();
    storageFloat[tag+"_vProdX"	    ] = event[idx].xProd();
    storageFloat[tag+"_vProdY"	    ] = event[idx].yProd();
    storageFloat[tag+"_vProdZ"	    ] = event[idx].zProd();
    storageFloat[tag+"_vDecayX"	] = event[idx].xDec();
    storageFloat[tag+"_vDecayY"	] = event[idx].yDec();
    storageFloat[tag+"_vDecayZ"	] = event[idx].zDec();
    
    Float_t dTemp = sqrt( event[idx].xProd() *  event[idx].xProd() + event[idx].yProd()* event[idx].yProd());
    storageFloat[tag+"_vProdLXY"    ] = dTemp;
    dTemp = sqrt( event[idx].xProd() *  event[idx].xProd() + event[idx].yProd()* event[idx].yProd() +  event[idx].zProd()* event[idx].zProd());
    storageFloat[tag+"_vProdR"      ] = dTemp; 
    
    dTemp = sqrt( event[idx].xDec() *  event[idx].xDec() + event[idx].yDec()* event[idx].yDec());
    storageFloat[tag+"_vDecayLXY"  ] = dTemp; 
    dTemp = sqrt( event[idx].xDec() *  event[idx].xDec() + event[idx].yDec()* event[idx].yDec() +  event[idx].zDec()* event[idx].zDec());
    storageFloat[tag+"_vDecayR"    ] = dTemp; 
    
    dTemp = sqrt(( event[idx].xDec() - event[idx].xProd() ) *  (event[idx].xDec() -event[idx].xProd()) + (event[idx].yDec()-event[idx].yProd())*( event[idx].yDec()-event[idx].yProd()));
    storageFloat[tag+"_decayDLXY"    ] = dTemp; 

    dTemp = sqrt( ( event[idx].xDec() - event[idx].xProd() ) *  (event[idx].xDec() -event[idx].xProd()) + (event[idx].yDec()-event[idx].yProd())*( event[idx].yDec()-event[idx].yProd()) + ( event[idx].zDec() - event[idx].zProd() )*(event[idx].zDec() -event[idx].zProd()) );
    storageFloat[tag+"_decayDR"     ] = dTemp; 
    storageFloat[tag+"_decayDZ"     ] =  abs(event[idx].zProd() - event[idx].zDec()); 
    
}
