package edu.rice.cs.bioinfo.programs.soranus;

import edu.rice.cs.bioinfo.library.programming.Proc;
import edu.rice.cs.bioinfo.library.programming.Proc1;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.soranus.controllers.AnalysisController;
import edu.rice.cs.bioinfo.programs.soranus.controllers.DataController;
import edu.rice.cs.bioinfo.programs.soranus.models.data.*;
import edu.rice.cs.bioinfo.programs.soranus.viewModels.WorkspaceVM;
import edu.rice.cs.bioinfo.programs.soranus.views.WorkspaceView;

import java.io.File;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/1/13
 * Time: 3:03 PM
 * To change this template use File | Settings | File Templates.
 */
abstract class ApplicationBase
{
    public void start()
    {
        SequencingsDataProvider<File, String, Exception> sequencingsDataProvider = makeSequencingsDataProvider();
        FirstPositiveDataProvider<File,String,Exception> firstPositiveDataProvider = makeFirstPositiveDataProvider();
        TraceDataProvider<File,String,Exception> traceDataProvider = makeTraceDataProvider();


        WorkspaceVM<DataController.FileDataRecord, AnalysisController.AnalysisRecord> vm =
                new WorkspaceVM<DataController.FileDataRecord, AnalysisController.AnalysisRecord>();
        DataController dc = new DataController(vm);
        AnalysisController<File,File,File,String,Exception> ac =
                new AnalysisController<File,File,File, String, Exception>(vm, sequencingsDataProvider, traceDataProvider,
                        firstPositiveDataProvider);
        WorkspaceView view = makeWorkspaceView(vm);
        configureViewForDataController(view, dc);
        configureViewForAnalysisController(view, ac);
        view.startView();

    }

    private TraceDataProvider<File, String, Exception> makeTraceDataProvider() {

        return new TraceDataProviderFile<String>() {
            @Override
            protected String makeSourceId(String sourceIdText) {
                return sourceIdText;  //To change body of implemented methods use File | Settings | File Templates.
            }
        };

    }

    private FirstPositiveDataProvider<File, String, Exception> makeFirstPositiveDataProvider() {

        return new FirstPositiveDataProviderFile<String>() {
            @Override
            protected String makeSourceId(String sourceIdText) {
                return sourceIdText;
            }
        };

    }

    private SequencingsDataProvider<File, String, Exception> makeSequencingsDataProvider()
    {
        return new SequencingsDataProviderFile<String>()
        {
            @Override
            protected String makeSourceId(String sourceIdText) {
                return sourceIdText;
            }
        };
    }

    protected void configureViewForDataController(WorkspaceView<DataController.FileDataRecord,
            AnalysisController.AnalysisRecord> view,
                                                  final DataController dc)
    {
        view.getAddDataRequested().addObserver(new Proc1<File>() {
            public void execute(File input) {
                try
                {
                    dc.addDataFile(input);
                }
                catch(Exception e)
                {
                    reportException(e);
                }
            }
        });

        view.getDataRecordSelected().addObserver(new Proc1<DataController.FileDataRecord>()
        {
            public void execute(DataController.FileDataRecord input) {
                dc.dataRecordSelected(input);
            }
        });

        view.addCreateSequencingsDataRequestListener(new Proc() {
            public void execute() {
                try
                {
                    dc.createSequencingDataFileFromOpenVaalOutputFiles();
                }
                catch (Exception e)
                {
                    reportException(e);
                }
            }
        });
    }

    protected void configureViewForAnalysisController(WorkspaceView<DataController.FileDataRecord,
            AnalysisController.AnalysisRecord> view,
                                                      final AnalysisController<File,File,File,String,Exception> ac)
    {
        view.getSnitkinTransMapAnalysisRequested().addObserver(new Proc3<DataController.FileDataRecord, DataController.FileDataRecord,
                DataController.FileDataRecord>() {
            public void execute(DataController.FileDataRecord sequencings, DataController.FileDataRecord traces,
                                DataController.FileDataRecord firstPositive) {

                try
                {
                    ac.performSnitkinTransMapAnalysis(sequencings.File, traces.File, firstPositive.File);
                }
                catch(Exception e)
                {
                    reportException(e);
                }
            }
        });

        view.getNeighborJoiningAnalysisRequested().addObserver(new Proc1<DataController.FileDataRecord>()
        {
            public void execute(DataController.FileDataRecord sequencings)
            {
                try {
                    ac.performNeighborJoiningAnalysis(sequencings.File);
                } catch (Exception e) {
                    reportException(e);
                }
            }
        });

        view.getDetectRecombAnalysisRequested().addObserver(new Proc1<DataController.FileDataRecord>()
        {
            public void execute(DataController.FileDataRecord sequencings)
            {
                try {
                    ac.performRecombDetectionUnderInfSite(sequencings.File);
                } catch (Exception e) {
                    reportException(e);
                }
            }
        });

        view.getMinSpanTreesSnpRequested().addObserver(new Proc1<DataController.FileDataRecord>()
        {
            public void execute(DataController.FileDataRecord sequencings)
            {
                try
                {
                    ac.performMinSpanTreeAnalysisSnp(sequencings.File);
                }
                catch(Exception e)
                {
                    reportException(e);
                }
            }
        });

        view.getMinSpanTreesMax2SnpRequested().addObserver(new Proc1<DataController.FileDataRecord>()
        {
            public void execute(DataController.FileDataRecord sequencings)
            {
                try
                {
                    ac.performMinSpanTreeAnalysisSnpMax2(sequencings.File);
                }
                catch(Exception e)
                {
                    reportException(e);
                }
            }
        });

        view.getAnalysisRecordSelected().addObserver(new Proc1<AnalysisController.AnalysisRecord>() {
            public void execute(AnalysisController.AnalysisRecord record) {
                ac.analysisRecordSelected(record);
            }
        });


    }

    private void reportException(Exception e) {
        throw new RuntimeException(e);
    }

    protected abstract WorkspaceView makeWorkspaceView(WorkspaceVM vm);


}
