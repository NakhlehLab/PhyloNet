package edu.rice.cs.bioinfo.programs.soranus;

import edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.Snitkin2012.SnitkinEdge;
import edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.TransMapResult;
import edu.rice.cs.bioinfo.library.programming.FuncEx;
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
import edu.rice.cs.bioinfo.programs.soranus.viewModels.*;
import edu.rice.cs.bioinfo.programs.soranus.views.WorkspaceView;
import org.apache.commons.io.FileUtils;
import org.joda.time.LocalDate;
import org.xml.sax.SAXException;

import javax.xml.parsers.ParserConfigurationException;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/6/13
 * Time: 4:00 PM
 * To change this template use File | Settings | File Templates.
 */
abstract class ControllerBase
{
    private edu.rice.cs.bioinfo.programs.soranus.viewModels.WorkspaceVM _workspaceViewModel;

    private WorkspaceView _workspaceView;

    private Map<String,File> _explorerTitleToDataFile = new HashMap<String, File>();

    private Map<edu.rice.cs.bioinfo.programs.soranus.viewModels.WorkspaceVM.AnalysisRecord,DocumentVM> _analysisRecordToAnalysisVM = new HashMap<WorkspaceVM.AnalysisRecord, DocumentVM>();

    private Map<WorkspaceVM.DataRecord,DocumentVM> _dataRecordToDataVM = new HashMap<WorkspaceVM.DataRecord, DocumentVM>();

    public void startView()
    {
        _workspaceViewModel = new WorkspaceVM();
        _workspaceView = makeWorkspaceView(_workspaceViewModel);
        _workspaceView.addDataAddRequestListener(new Proc1<File>() {
            public void execute(File dataFile) {
                try
                {
                    addDataFile(dataFile);
                }
                catch (Exception e)
                {
                    reportException(e);
                }
            }
        });
        _workspaceView.addSnitkinTransMapAnalysisRequestedListener(new Proc3<String, String, String>() {
            public void execute(final String sequencingsTitle, final String traceTitle, final String firstPositiveTitle) {

                try
                {
                    performAnalysis(new FuncEx<DocumentVM, Exception>() {
                        public DocumentVM execute() throws Exception {
                            return performSnitkinTransMapAnalysis(sequencingsTitle, traceTitle, firstPositiveTitle);
                        }
                    }, "Snitkin Trans Map");
                }
                catch (Exception e)
                {
                    reportException(e);
                }
            }
        });
        _workspaceView.addNeighborJoiningAnalysisRequestedListener(new Proc1<String>() {
            public void execute(final String sequencingsTitle) {
                try
                {
                    performAnalysis(new FuncEx<DocumentVM,Exception>() {
                        public DocumentVM execute() throws Exception{
                            return  performNeighborJoiningAnalysis(sequencingsTitle);
                        }
                    }, "Neighbor Joining");

                }
                catch (Exception e)
                {
                    reportException(e);
                }
            }
        });
        _workspaceView.addAnalysisRecordSelectedListener(new Proc1<WorkspaceVM.AnalysisRecord>() {
            public void execute(WorkspaceVM.AnalysisRecord record) {
                DocumentVM analysisVM = _analysisRecordToAnalysisVM.get(record);
                _workspaceViewModel.setFocusDocument(analysisVM);
            }
        });
        _workspaceView.addDataRecordSelectedListener(new Proc1<WorkspaceVM.DataRecord>() {
            public void execute(WorkspaceVM.DataRecord record) {
                DocumentVM dataVM = _dataRecordToDataVM.get(record);
                _workspaceViewModel.setFocusDocument(dataVM);
            }
        });



        _workspaceView.startView();
    }

    protected void reportException(Exception e)
    {
        int j = 0;
    }


    protected abstract WorkspaceView makeWorkspaceView(WorkspaceVM vm);

    private void addDataFile(final File dataFile) throws IOException
    {

        KnownDatafileFormat format = new DatafileRecogniser().recognise(dataFile);

        String contents = FileUtils.readFileToString(dataFile);
        XMLDataVM vm = new XMLDataVM(contents);

        WorkspaceVM.DataRecord record = format.execute(new KnownDatafileFormatAlgo<WorkspaceVM.DataRecord, Void, RuntimeException>() {
            public WorkspaceVM.DataRecord forSequencings(KnownDatafileFormat format, Void input) throws RuntimeException {
                String explorerTitle = dataFile.getName();
                _explorerTitleToDataFile.put(explorerTitle, dataFile);
                return _workspaceViewModel.addSequencingsData(explorerTitle);
            }

            public WorkspaceVM.DataRecord forTraces(KnownDatafileFormat format, Void input) throws RuntimeException {
                String explorerTitle = dataFile.getName();
                _explorerTitleToDataFile.put(explorerTitle, dataFile);
                return _workspaceViewModel.addTraceData(explorerTitle);
            }

            public WorkspaceVM.DataRecord forFirstPositive(KnownDatafileFormat format, Void input) throws RuntimeException {
                String explorerTitle = dataFile.getName();
                _explorerTitleToDataFile.put(explorerTitle, dataFile);
                return _workspaceViewModel.addFirstPositiveData(explorerTitle);
            }
        }, null);

        _dataRecordToDataVM.put(record, vm);
        _workspaceViewModel.setFocusDocument(vm);
    }

    private void performAnalysis(FuncEx<DocumentVM,Exception> analysisCommand, String title) throws Exception
    {
        DocumentVM analysisVM = analysisCommand.execute();
        WorkspaceVM.AnalysisRecord record = _workspaceViewModel.addAnalysis(title);
        _analysisRecordToAnalysisVM.put(record, analysisVM);
        _workspaceViewModel.setFocusDocument(analysisVM);
    }

    private DocumentVM performSnitkinTransMapAnalysis(String sequencingsTitle, String traceTitle, String firstPositiveTitle) throws IOException, SAXException, ParserConfigurationException {

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

            return transMapViewModel;


    }

    private DocumentVM performNeighborJoiningAnalysis(String sequencingsTitle) throws IOException, SAXException, ParserConfigurationException {


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


            NeighborJoiningVM<Sequencing,NeighborJoining.Edge> joinVM = new NeighborJoiningVM<Sequencing,NeighborJoining.Edge>(joinTree.Nodes, joinTree.Edges)
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

                @Override
                public String getEdgeLabel(NeighborJoining.Edge edge) {
                    return edge.Distance.toString();
                }
            };
            return joinVM;


    }
}
