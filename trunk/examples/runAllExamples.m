function runAllExamples()

cd(['hessians',filesep,'logsumexp']);
main;
cd(['..',filesep,'..']);
drawnow ; input ('hit enter to continue:') ;


cd(['jachesvecprods']);
main;
cd(['..']);
drawnow ; input ('hit enter to continue:') ;


cd(['jacobians',filesep,'arrowhead']);
main;
cd(['..',filesep,'..']);
drawnow ; input ('hit enter to continue:') ;


cd(['jacobians',filesep,'polydatafit']);
main;
cd(['..',filesep,'..']);
drawnow ; input ('hit enter to continue:') ;


cd(['optimization',filesep,'fminconEx']);
main;
cd(['..',filesep,'..']);
drawnow ; input ('hit enter to continue:') ;


cd(['optimization',filesep,'fminuncEx']);
main;
cd(['..',filesep,'..']);
drawnow ; input ('hit enter to continue:') ;


cd(['optimization',filesep,'fsolveEx']);
main;
cd(['..',filesep,'..']);
drawnow ; input ('hit enter to continue:') ;


cd(['optimization',filesep,'ipoptEx']);
gl2main;
cd(['..',filesep,'..']);
drawnow ; input ('hit enter to continue:') ;

cd(['optimization',filesep,'fsolveEx']);
main;
cd(['..',filesep,'..']);
drawnow ; input ('hit enter to continue:') ;


cd(['optimization',filesep,'vectorized',filesep,'brachistochrone']);
main_basic_1stderivs;
drawnow ; input ('hit enter to continue:') ;
main_basic_2ndderivs;
drawnow ; input ('hit enter to continue:') ;
main_vect_1stderivs;
drawnow ; input ('hit enter to continue:') ;
main_vect_2ndderivs;
cd(['..',filesep,'..',filesep,'..']);
drawnow ; input ('hit enter to continue:') ;
close all;


cd(['optimization',filesep,'vectorized',filesep,'minimumclimb']);
main_1stderivs_nonvect;
drawnow ; input ('hit enter to continue:') ;
main_1stderivs_vect;
drawnow ; input ('hit enter to continue:') ;
main_2ndderivs_nonvect;
drawnow ; input ('hit enter to continue:') ;
main_2ndderivs_vect;
cd(['..',filesep,'..',filesep,'..']);
drawnow ; input ('hit enter to continue:') ;
close all;


cd(['stiffodes',filesep,'brusselator']);
main;
cd(['..',filesep,'..']);
drawnow ; input ('hit enter to continue:') ;


cd(['stiffodes',filesep,'burgers']);
main;
cd(['..',filesep,'..']);
drawnow ; input ('hit enter to continue:') ;


cd(['stiffodes',filesep,'DCALcontrol']);
main;
cd(['..',filesep,'..']);
drawnow ; input ('hit enter to continue:') ;

close all;
fprintf(1,'Successfully ran all adigator example problems.\n');
end
