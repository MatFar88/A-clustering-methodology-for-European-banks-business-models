cd('C:\Users\Calif\OneDrive\Documenti\Passato\Desktop\Backup PC\Documenti\Expansion Drive\Documents\Documents\Università\Post-Doc\spec');
%%
z=mvnrnd(repmat(0,100,1),diag(repmat(1,100,1)),100)

%%
tfkm_try=tfkm(z,0.1,2:3,2:3,0)
tfkm_try_2=tfkm(z,0.1,2:3,4,0)
tfkm_try_3=tfkm(z,0.1,4,2:3,0)
tfkm_try_4=tfkm(z,0.1,4,4,0)
tfkm_try_5=tfkm(z,0.1,2,4,0)
tfkm_try_6=tfkm(z,0.1,2:3,2:4,0,0,100)
tfkm_try_6{1,1}
tfkm_try_6{1,2}
tfkm_try_6{1,15}
tfkm_try_6{1,22}

%%

%tfkm_try_alpha=tfkm_alpha(z,[0 0.1 0.2],2:3,2:3,0)
tfkm_try_alpha=tfkm_alpha(z,0.1,2:3,2:3,0)
tfkm_try_alpha=tfkm_alpha(z,0,2:3,2:3,0)
tfkm_try_alpha=tfkm_alpha(z,[0 0.1 0.2],2:3,2:3,0)
tfkm_try_alpha_2=tfkm_alpha(z,0.1,2:3,4,0)
tfkm_try_alpha_2=tfkm_alpha(z,[0 0.1 0.2],2:3,4,0)
tfkm_try_alpha_2=tfkm_alpha(z,[0],2:3,4,0)
tfkm_try_alpha_3=tfkm_alpha(z,[0 0.1 0.2],4,2:3,0)
tfkm_try_alpha_3=tfkm_alpha(z,0.1,4,2:3,0)
tfkm_try_alpha_3=tfkm_alpha(z,[0],4,2:3,0)
tfkm_try_alpha_4=tfkm_alpha(z,0.1,4,4,0)
tfkm_try_alpha_4=tfkm_alpha(z,[0 0.1 0.2],4,4,0)
tfkm_try_alpha_4=tfkm_alpha(z,[0],4,4,0)
tfkm_try_alpha_5=tfkm_alpha(z,0.1,2,4,0)
tfkm_try_alpha_5=tfkm_alpha(z,[0 0.1 0.2],2,4,0)
tfkm_try_alpha_5=tfkm_alpha(z,[0],2,4,0)
tfkm_try_alpha_6_1=tfkm_alpha(z,0.1,2:3,2:4,0)
tfkm_try_alpha_6_2=tfkm_alpha(z,0,2:3,2:4,0)
tfkm_try_alpha_6_3=tfkm_alpha(z,[0 0.1 0.2],2:3,2:4,0,0,100)