package edu.rice.cs.bioinfo.programs.soranus;

import edu.rice.cs.bioinfo.library.programming.Proc1;
import edu.rice.cs.bioinfo.library.programming.Proc2;
import edu.rice.cs.bioinfo.library.programming.Proc3;
import edu.rice.cs.bioinfo.programs.soranus.models.analysis.Sequencing;
import edu.rice.cs.bioinfo.programs.soranus.models.factories.FirstPositiveMapFactory;
import edu.rice.cs.bioinfo.programs.soranus.models.factories.PatientTraceMapFactory;
import edu.rice.cs.bioinfo.programs.soranus.models.factories.SequencingsMapFactory;
import edu.rice.cs.bioinfo.programs.soranus.models.fileRecogniser.DatafileRecogniser;
import edu.rice.cs.bioinfo.programs.soranus.models.fileRecogniser.KnownDatafileFormat;
import edu.rice.cs.bioinfo.programs.soranus.models.fileRecogniser.KnownDatafileFormatAlgo;
import edu.rice.cs.bioinfo.programs.soranus.viewModels.WorkspaceVM;
import edu.rice.cs.bioinfo.programs.soranus.views.swing.Workspace;
import org.apache.commons.io.FileUtils;
import org.joda.time.LocalDate;

import javax.swing.*;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

class DefaultController
{
    private WorkspaceVM _workspaceViewModel;

    private Workspace _workspaceView;

    private Map<String,File> _explorerTitleToDataFile = new HashMap<String, File>();

    DefaultController()
    {
        SwingUtilities.invokeLater(new Runnable()
        {
            public void run()
            {

                _workspaceViewModel = new WorkspaceVM();
                _workspaceView = new Workspace(_workspaceViewModel);
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
                    public void execute(String sequencingsTitle, String traceTitle, String firstPositiveTitle)
                    {
                        doSnitkinTransMapAnalysisRequested(sequencingsTitle, traceTitle, firstPositiveTitle);

                    }
                });

            }
        });
    }

    private void doSnitkinTransMapAnalysisRequested(String sequencingsTitle, String traceTitle, String firstPositiveTitle)
    {
        try
        {
            File sequencingsFile = _explorerTitleToDataFile.get(sequencingsTitle);
            File traceFile = _explorerTitleToDataFile.get(traceTitle);
            File firstPositiveFile = _explorerTitleToDataFile.get(firstPositiveTitle);

            Map<Sequencing,Integer> sequeincingToPatient = new SequencingsMapFactory<Integer>()
            {

                @Override
                protected Integer makeId(String sourceIdText) {
                    return new Integer(sourceIdText);
                }
            }.makeMapFromXMLFile(sequencingsFile);

            Map<Integer,LocalDate> patientToFirstPositiveDate = new FirstPositiveMapFactory<Integer,LocalDate>()
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

            Map<Integer, Map<LocalDate, Integer>> patientTraces = new PatientTraceMapFactory<Integer,LocalDate,Integer>()
            {
                @Override
                protected Integer makeLocation(String locationText) {
                    return new Integer(locationText);  //To change body of implemented methods use File | Settings | File Templates.
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

            int j = 0;
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
