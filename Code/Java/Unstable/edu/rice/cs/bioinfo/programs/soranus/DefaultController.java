package edu.rice.cs.bioinfo.programs.soranus;

import edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.Snitkin2012.NIHOutbreakDataTestBase;
import edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.Snitkin2012.SnitkinEdge;
import edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.Snitkin2012.SnitkinTransMapInferrerTemplateDoubleDistance;
import edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.TransMapResult;
import edu.rice.cs.bioinfo.library.programming.Proc1;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.library.programming.Tuple;
import edu.rice.cs.bioinfo.programs.soranus.models.analysis.NeighborJoining;
import edu.rice.cs.bioinfo.programs.soranus.models.analysis.Sequencing;
import edu.rice.cs.bioinfo.programs.soranus.models.analysis.SnitkinTransMapInferrer;
import edu.rice.cs.bioinfo.programs.soranus.models.factories.FirstPositiveMapFactory;
import edu.rice.cs.bioinfo.programs.soranus.models.factories.PatientTraceMapFactory;
import edu.rice.cs.bioinfo.programs.soranus.models.factories.SequencingsMapFactory;
import edu.rice.cs.bioinfo.programs.soranus.models.fileRecogniser.DatafileRecogniser;
import edu.rice.cs.bioinfo.programs.soranus.models.fileRecogniser.KnownDatafileFormat;
import edu.rice.cs.bioinfo.programs.soranus.models.fileRecogniser.KnownDatafileFormatAlgo;
import edu.rice.cs.bioinfo.programs.soranus.viewModels.NeighborJoiningVM;
import edu.rice.cs.bioinfo.programs.soranus.viewModels.TransMapVM;
import edu.rice.cs.bioinfo.programs.soranus.viewModels.WorkspaceVM;
import edu.rice.cs.bioinfo.programs.soranus.views.swing.NeighborJoiningView;
import edu.rice.cs.bioinfo.programs.soranus.views.swing.WorkspaceView;
import org.joda.time.LocalDate;

import javax.swing.*;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

class DefaultController
{
    private WorkspaceVM _workspaceViewModel;

    private WorkspaceView _workspaceView;

    private Map<String,File> _explorerTitleToDataFile = new HashMap<String, File>();

    DefaultController()
    {
        SwingUtilities.invokeLater(new Runnable()
        {
            public void run()
            {

                _workspaceViewModel = new WorkspaceVM();
                _workspaceView = new WorkspaceView(_workspaceViewModel);
                _workspaceView.addAddDataListener(new Proc1<File>() {
                    public void execute(File dataFile)
                    {
                        try
                        {
                            onAddData(dataFile);
                        }
                        catch(IOException e)
                        {
                            reportException(e) ;
                        }

                    }
                });
                _workspaceView.addSnitkinTransMapAnalysisRequestedListener(new Proc3<String, String, String>() {
                    public void execute(String sequencingsTitle, String traceTitle, String firstPositiveTitle) {
                        performSnitkinTransMapAnalysis(sequencingsTitle, traceTitle, firstPositiveTitle);

                    }
                });
                _workspaceView.addNeighborJoiningAnalysisRequestedListener(new Proc1<String>() {
                    public void execute(String sequencingsTitle) {
                        performNeighborJoiningAnalysis(sequencingsTitle);
                    }
                });

            }
        });
    }

    private void performNeighborJoiningAnalysis(String sequencingsTitle) {

        try
        {
            File sequencingsFile = _explorerTitleToDataFile.get(sequencingsTitle);

            final Map<Sequencing,Integer> sequeincingToPatient = new SequencingsMapFactory<Integer>()
            {

                @Override
                protected Integer makeId(String sourceIdText) {
                    return new Integer(sourceIdText);
                }
            }.makeMapFromXMLFile(sequencingsFile);

            NeighborJoining nj = new NeighborJoining();
            NeighborJoining.Graph joinTree = nj.performJoin(sequeincingToPatient.keySet());
            NeighborJoiningVM<Sequencing,NeighborJoining.Edge> joinVM = new NeighborJoiningVM<Sequencing,NeighborJoining.Edge>(joinTree.Edges)
            {
                @Override
                public Tuple<Sequencing, Sequencing> getNodesOfEdge(NeighborJoining.Edge edge) {
                    return new Tuple<Sequencing, Sequencing>(edge.Node1, edge.Node2);
                }

                @Override
                public String getNodeLabel(Sequencing node) {
                    if(sequeincingToPatient.containsKey(node))
                    {
                        return sequeincingToPatient.get(node).toString();
                    }
                    else
                    {
                        return "";
                    }
                }
            };
            joinVM.setRoot(nj.getLastNodeCreated());
            _workspaceViewModel.setFocusDocument(joinVM);
        }
        catch(Exception e)
        {
            reportException(e);
        }
    }

    private void performSnitkinTransMapAnalysis(String sequencingsTitle, String traceTitle, String firstPositiveTitle)
    {
        try
        {
            File sequencingsFile = _explorerTitleToDataFile.get(sequencingsTitle);
            File traceFile = _explorerTitleToDataFile.get(traceTitle);
            File firstPositiveFile = _explorerTitleToDataFile.get(firstPositiveTitle);

            final Map<Sequencing,Integer> sequeincingToPatientOuter = new SequencingsMapFactory<Integer>()
            {

                @Override
                protected Integer makeId(String sourceIdText) {
                    return new Integer(sourceIdText);
                }
            }.makeMapFromXMLFile(sequencingsFile);

            final Map<Integer,LocalDate> patientToFirstPositiveDateOuter = new FirstPositiveMapFactory<Integer,LocalDate>()
            {
                @Override
                protected LocalDate makeDate(String dateText) {
                    return new LocalDate(dateText);
                }

                @Override
                protected Integer makeId(String sourceIdText) {
                    return new Integer(sourceIdText);
                }

            }.makeMapFromXMLFile(firstPositiveFile);

            final Map<Integer, Map<LocalDate, Object>> patientTracesOuter = new PatientTraceMapFactory<Integer,LocalDate,Object>()
            {
                @Override
                protected Object makeLocation(String locationText) {
                    return locationText;  //To change body of implemented methods use File | Settings | File Templates.
                }

                @Override
                protected LocalDate makeDate(String dateText) {
                    return new LocalDate(dateText);  //To change body of implemented methods use File | Settings | File Templates.
                }

                @Override
                protected Integer makeId(String sourceIdText) {
                    return new Integer(sourceIdText);  //To change body of implemented methods use File | Settings | File Templates.
                }
            }.makeMapFromXMLFile(traceFile);

            final TransMapResult<SnitkinEdge<Integer,Double>> transResult = new SnitkinTransMapInferrer(patientTracesOuter, patientToFirstPositiveDateOuter, sequeincingToPatientOuter).inferMaps(1);
            final Set<SnitkinEdge<Integer,Double>> transMap = transResult.getSolutions().iterator().next();

            TransMapVM<Integer, SnitkinEdge<Integer,Double>> transMapViewModel = new TransMapVM<Integer, SnitkinEdge<Integer, Double>>(transMap) {
                @Override
                public Integer getSource(SnitkinEdge<Integer, Double> edge) {
                    return edge.getSource();
                }

                @Override
                public Integer getDestination(SnitkinEdge<Integer, Double> edge) {
                    return edge.getDestination();
                }

                @Override
                public String getNodeLabel(Integer integer) {
                    return integer.toString();
                }
            };

            _workspaceViewModel.setFocusDocument(transMapViewModel);
        }
        catch (Exception e)
        {
            reportException(e);
        }
    }

    private void reportException(Exception e)
    {

    }

    private void onAddData(final File dataFile) throws IOException
    {

        KnownDatafileFormat format = new DatafileRecogniser().recognise(dataFile);
        format.execute(new KnownDatafileFormatAlgo<Void, Void, RuntimeException>() {
            public Void forSequencings(KnownDatafileFormat format, Void input) throws RuntimeException {
                String explorerTitle = dataFile.getName();
                _explorerTitleToDataFile.put(explorerTitle, dataFile);
                _workspaceViewModel.addSequencingsData(explorerTitle);
                return null;
            }

            public Void forTraces(KnownDatafileFormat format, Void input) throws RuntimeException {
                String explorerTitle = dataFile.getName();
                _explorerTitleToDataFile.put(explorerTitle, dataFile);
                _workspaceViewModel.addTraceData(explorerTitle);
                return null;
            }

            public Void forFirstPositive(KnownDatafileFormat format, Void input) throws RuntimeException {
                String explorerTitle = dataFile.getName();
                _explorerTitleToDataFile.put(explorerTitle, dataFile);
                _workspaceViewModel.addFirstPositiveData(explorerTitle);
                return null;
            }
        }, null);
    }

    public void startView()
    {
        SwingUtilities.invokeLater(new Runnable() {
            public void run() {

                _workspaceView.setExtendedState(_workspaceView.getExtendedState() | JFrame.MAXIMIZED_BOTH);
                _workspaceView.setVisible(true);
            }
        });
    }

}
