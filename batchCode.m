clear
%%%%%
% add everything to path if any errors
%%%%%
cd("/Users/krithika/Library/CloudStorage/Box-Box/kk_june/2023_codes_kk/")

% Add programs folders to path
dynamo = "/Users/krithika/Library/CloudStorage/Box-Box/kk_june/2023_codes_kk/march/dynamo-master"
dynamo_libsbml = "/Users/krithika/Library/CloudStorage/Box-Box/kk_june/2023_codes_kk/march/dynamo-master/codes_libSBML"
lib_sbml = ("march/libSBML-5.19.0-matlab-binaries-2/")

% setting paths and directories to extract BioModels and save results
sbml_toolbox = "march/SBMLToolbox-4.1.0/"

% BioModels = "/Users/krithika/Library/CloudStorage/Box-Box/kk_june/2023_codes_kk/march/BioModels/"
BioModels = 'march/BioModels/'

% destination folder to save results
res.odes = "/Users/krithika/Library/CloudStorage/Box-Box/kk_june/2023_codes_kk/march/res.ODEs" 
res.ftabs = "/Users/krithika/Library/CloudStorage/Box-Box/kk_june/2023_codes_kk/march/res.ftabs"
res.jtabs = "/Users/krithika/Library/CloudStorage/Box-Box/kk_june/2023_codes_kk/march/res.jtabs"
res.fhattabs = "/Users/krithika/Library/CloudStorage/Box-Box/kk_june/2023_codes_kk/march/res.fhattabs"

res.maxSpeciesPerturbed = "/Users/krithika/Library/CloudStorage/Box-Box/kk_june/2023_codes_kk/march/res.maxSpeciesPerturbed"
res.mds = "/Users/krithika/Library/CloudStorage/Box-Box/kk_june/2023_codes_kk/march/mds"

plots.sims = "/Users/krithika/Library/CloudStorage/Box-Box/kk_june/2023_codes_kk/march/plots.sims"
plots.fValues = "/Users/krithika/Library/CloudStorage/Box-Box/kk_june/2023_codes_kk/march/plots.fValues"
plots.jacs = "/Users/krithika/Library/CloudStorage/Box-Box/kk_june/2023_codes_kk/march/plots.jacs"


addpath(dynamo,dynamo_libsbml,lib_sbml,BioModels,sbml_toolbox,res.odes,res.ftabs,res.jtabs,res.maxSpeciesPerturbed,plots.sims,plots.fValues,plots.jacs)
addpath(res.fhattabs)
set(0,'DefaultFigureVisible','off')

% grab 87 biomodels 

% Get a list of all files in the folder with the desired file name pattern.
filePattern = fullfile(BioModels, '*.xml'); % Change to whatever pattern you need.
theFiles = dir(filePattern);

% k= # 1-1,3-13,5,6-15,7,10,11(good),12-38,13-046,14-061,15-069,16-070,17-85,19-86,65-398,86-600,87-602 (good)
% errors with k= # 2(biomod9) -> fix manually line 230 , 4, 6, 8, 9 (mostly issues with tolerance)
% modError = [6,8,9];

%  broke at k= 4 - > Expected DATA to be real.


% 12,15 20, 30, 61, 66, 77, 85
% [1,3,5,7,10,11,65]
for k = 1 : length(theFiles)
       try     
            k=1
            baseFileName = theFiles(k).name;
            fullFileName = fullfile(theFiles(k).folder, baseFileName);
            [p,f,e]=fileparts(fullFileName);
            filename=fullfile(f);
            fprintf(1, 'BioModel %s\n', baseFileName);
              
            
            model = baseFileName;
            model_id = filename;
            fprintf(model);
            cd("/Users/krithika/Library/CloudStorage/Box-Box/kk_june/2023_codes_kk/march/BioModels/");
            [SBMLModel, errors] = TranslateSBML((model),1,0);
            nSpecies = length(SBMLModel.species);
            
            
            Name = '';
            if (SBMLModel.SBML_level == 1)
                Name = SBMLModel.name;
            else
                if (isempty(SBMLModel.id))
                    Name = SBMLModel.name;
                else
                    Name = SBMLModel.id;
                end
            end
                
            if (length(Name) > 63)
                Name = Name(1:60);
            end
            
            model_name = Name;
            model_handle = [ '@' model_name ];
            
                       
            SpeciesNames = GetSpecies(SBMLModel);
            
            % uncomment this out below if you want fullnames
            
%             if (SBMLModel.SBML_level == 1 || SBMLModel.SBML_level == 2) 
%             SpeciesFullNames = {}
%                 for x = 1:nSpecies
%                     temp = SBMLModel.species(x).name;
%                     SpeciesFullNames =  [SpeciesFullNames temp];
%                 end   
%             end
%             
%           
%             
%             if (isempty(SpeciesFullNames) == 0 && k~=9)
%                 SpeciesNames = SpeciesFullNames;
%             end    
%                 
            
            [VarParams, VarInitValues] = GetVaryingParameters(SBMLModel); % empty
            [ParNames, Parvalues] = GetAllParametersUnique(SBMLModel);
            VarNames = [ SpeciesNames VarParams ]; % species names
            
            cd(res.ods);
            WriteODEFunction(SBMLModel);
            
            
            %% Steady-state and Jacobian
            % initial concentrations
            x0 = eval(model_name);
            
            
            initSpeciesData = array2table(x0(1:nSpecies));
            initSpeciesData.Properties.RowNames = VarNames(1:nSpecies);
            
            speciesData = addvars(initSpeciesData);
            speciesData.Properties.VariableNames{1} = 'initialData';
            
            [t1,x1] = ode23tb(eval(model_handle), [0, 200], x0);
            xss1 = x1(end,:);% + randn() * 1e-15;
            F0 = array2table(max(x1(end,:),0));         
            F0 = renamevars(F0,1:width(F0),VarNames);
            
            
            % add simulated data to table
            speciesData = addvars(speciesData,transpose(xss1(:,1:nSpecies)));
            speciesData.Properties.VariableNames{2} = 'simulatedData (F0)';
            
            % plot simulation 1 metab change
            x1Plot = plot(t1, log(x1(:, 1:length(SpeciesNames(1:nSpecies)))));
%             ylim([min(min(x1)) max(max(x1))]);
            lgd = legend(SpeciesNames(1:nSpecies));
            xlabel('Time')
            ylabel('log(Species)');
            title(['BioModel: ' model_id(13:15)]);
            lgd.FontSize = 8;
            set(gca,"FontSize",12);
            
            cd(res.mds);
            simplot1 = ['Simulation1:' model_id(13:15) '.eps'];
            set(gcf, 'PaperUnits', 'inches');
            x_width=15 ;y_width=7.5;
            set(gcf, 'PaperPosition', [0 0 x_width y_width]); 
            saveas(gcf,simplot1,'epsc');


            % loop to adjust initial concentration as per simulation back into the
            % model 
            for j = 1:length(xss1)  
%                disp(SBMLModel.species(i).initialAmount) 
                 SBMLModel.species(j).initialConcentration = xss1(j); 
            end
             
            
          
            %%
%             %perturbation of the highest concentration species 
                temp = []
                for x = 1:nSpecies
                    temp = [temp SBMLModel.species(x).initialConcentration];
                end
                [maxSpecValue,maxSpecIndex] = max(temp(:))
   
            % grab species value to perturb # value in species(I) to set to 0 
            maxSpecName = SBMLModel.species(maxSpecIndex).id;
            maxSpecName = (maxSpecName);
            
%             maxSpecName = reshape([maxSpecName],[],3);
            
            cd(res.mds)
            
            M = table(maxSpecIndex, string(maxSpecName),maxSpecValue, 'VariableNames', { 'Index', 'Name','Value'} );
            writetable(M,['maxSpeciesPerturbed',model_id(13:15),'.txt'], 'WriteRowNames', true);
            
           
            
             %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            cd(res.mds)
 
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            %% simple perturbation
            
            
            % perturb #  value in species(I) set to 0 
            SBMLModel.species(maxSpecIndex).initialConcentration = 0   
            x0_perturbed = [SBMLModel.species.initialConcentration];

            
            
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            % commmon step for getting delta_hat and x0_perturbed back into
            % the SBML model
            x0_perturbed = transpose(x0_perturbed);
       
            
            [t2,x2] = ode23tb(eval(model_handle), [0, 200],x0_perturbed);
            xss2 = x2(end,:);

            Ft = array2table(max(x2(end,:),0));           
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
            
            speciesData = addvars(speciesData,x0_perturbed(1:nSpecies,:));
            speciesData.Properties.VariableNames{3} = 'perturbedData';
           
            speciesData = addvars(speciesData,transpose(xss2(:,1:nSpecies)));
            speciesData.Properties.VariableNames{4} = 'simulatedData (Ft)';
                       
            cd(res.mds);
            % plot simulation 2 metab change
%             x2Plot = semilogy(t2, x2(:, 1:length(SpeciesNames(1:nSpecies))));    
            x2Plot = plot(t2, log(x2(:, 1:length(SpeciesNames(1:nSpecies)))));    
%             ylim([min(min(x1)) max(max(x1))]);
            lgd = legend(SpeciesNames(1:nSpecies));
            xlabel('Time')
            ylabel('log(Species)');
            title(['BioModel: ' model_id(13:15)]);
            lgd.FontSize = 8;
            set(gca,"FontSize",12);
            
            simplot2 = ['Simulation2:' model_id(13:15) '.eps'];
            set(gcf, 'PaperUnits', 'inches');
            x_width=15 ;y_width=7.5;
            set(gcf, 'PaperPosition', [0 0 x_width y_width]); 
            saveas(gcf,simplot2,'epsc');
            
            
            
            
            F0 = rows2vars(F0(:,1:nSpecies));
            Ft = rows2vars(Ft(:,1:nSpecies));
%             fTable = addvars(F0,Ft);
            F0.Properties.VariableNames{2} = 'F0';
            Ft.Properties.VariableNames{2} = 'Ft';
            fTable = [F0(:,2) Ft(:,2)];
            fTable.Properties.RowNames = transpose(SpeciesNames);
                        
          
            
            cd(res.mds);
            writetable(fTable,['fTable',model_id(13:15),'.csv'],'WriteRowNames',true);
            
             
            J = compJacobian(@(x1) eval([ model_name '(10, x1)']),xss2);
            Jtab = array2table(J(1:nSpecies,1:nSpecies),'VariableNames', SpeciesNames, 'RowNames', SpeciesNames);
            Jtab_species = table2array(Jtab);
       
            cd(res.mds);
            writetable(Jtab, ['jac',model_id(13:15),'.csv'], 'WriteRowNames', true);
            
            cd(res.mds);
            J_map = HeatMap(Jtab_species);
            addXLabel(J_map,'Species','FontSize',12);
            addYLabel(J_map,'Species','FontSize',12);
            fig_heat_perturbed = figure;
            heat_p = plot(J_map,fig_heat_perturbed);
            heat_p.Title.String = 'Species Data Perturbed';
            set(gca,'XTickLabelRotation',45);
            heat_p.XTickLabels = SpeciesNames; 
            heat_p.YTickLabels = SpeciesNames; 
            saveas(gcf, [model_id(13:15),'_perturbedJac.eps'],'epsc');
            
            close
            
            cd(res.mds);
            newXlabels = [speciesData.Properties.RowNames]; 
            stackFig1 = stackedplot(speciesData,"Title","Species Data");
            ax = findobj(stackFig1.NodeChildren, 'Type','Axes');
            set([ax.YLabel],'Rotation',90,'HorizontalAlignment','Center', 'VerticalAlignment', 'Bottom')
            stackFig1.XLabel = "Species";
            set(ax,'XTick',(1:nSpecies),'XTickLabel',newXlabels);
            xtickangle(ax,45);
            set(stackFig1, 'FontSize', 14);
            set(gcf, 'PaperUnits', 'inches');
            x_width=15 ;y_width=11;
            set(gcf, 'PaperPosition', [0 0 x_width y_width]); 
            saveas(gcf, [model_id(13:15),'_Fvalues.eps'],'epsc');
   
               
        cd('../')
    catch
        fprintf('loop number %d failed\n',k);
         
    end
end     