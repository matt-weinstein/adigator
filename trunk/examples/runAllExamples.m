function runAllExamples()

cd(['hessians',filesep,'logsumexp']);
main;
cd(['..',filesep,'..']);


cd(['jachesvecprods']);
main;
cd(['..']);


cd(['jacobians',filesep,'arrowhead']);
main;
cd(['..',filesep,'..']);


cd(['jacobians',filesep,'polydatafit']);
main;
cd(['..',filesep,'..']);


cd(['optimization',filesep,'fminconEx']);
main;
cd(['..',filesep,'..']);


cd(['optimization',filesep,'fminuncEx']);
main;
cd(['..',filesep,'..']);


cd(['optimization',filesep,'fsolveEx']);
main;
cd(['..',filesep,'..']);


cd(['optimization',filesep,'ipoptEx']);
gl2main;
cd(['..',filesep,'..']);


cd(['optimization',filesep,'fsolveEx']);
main;
cd(['..',filesep,'..']);


cd(['optimization',filesep,'vectorized',filesep,'brachistochrone']);
main_basic_1stderivs;
main_basic_2ndderivs;
main_vect_1stderivs;
main_vect_2ndderivs;
cd(['..',filesep,'..',filesep,'..']);
close all;


cd(['optimization',filesep,'vectorized',filesep,'minimumclimb']);
main_1stderivs_nonvect;
main_1stderivs_vect;
main_2ndderivs_nonvect;
main_2ndderivs_vect;
cd(['..',filesep,'..',filesep,'..']);
close all;


cd(['stiffodes',filesep,'brusselator']);
main;
cd(['..',filesep,'..']);


cd(['stiffodes',filesep,'burgers']);
main;
cd(['..',filesep,'..']);


cd(['stiffodes',filesep,'DCALcontrol']);
main;
cd(['..',filesep,'..']);

close all;
clc;
fprintf(1,'Successfully ran all adigator example problems.\n');
end
