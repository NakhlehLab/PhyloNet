package edu.rice.cs.bioinfo.programs.soranus.controllers;

import edu.rice.cs.bioinfo.programs.soranus.models.factories.SequencingsXmlFromVaalOutFactory;
import edu.rice.cs.bioinfo.programs.soranus.models.fileRecogniser.DatafileRecogniser;
import edu.rice.cs.bioinfo.programs.soranus.models.fileRecogniser.KnownDatafileFormat;
import edu.rice.cs.bioinfo.programs.soranus.models.fileRecogniser.KnownDatafileFormatAlgo;
import edu.rice.cs.bioinfo.programs.soranus.viewModels.*;
import org.apache.commons.io.FileUtils;

import java.io.File;
import java.io.IOException;
import java.text.ParseException;
import java.util.HashMap;
import java.util.HashSet;
import java.util.Map;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 4/1/13
 * Time: 2:18 PM
 * To change this template use File | Settings | File Templates.
 */
public class DataController extends ControllerBase
{
    public class FileDataRecord extends DataRecord
    {
        public final File File;

        public final DocumentVM DocumentVM;

        protected FileDataRecord(String title, File file, DocumentVM documentVM) {
            super(title);
            File = file;
            DocumentVM = documentVM;
        }
    }

    private WorkspaceVM<FileDataRecord,?> _workspaceViewModel;

    private Set<File> _openVaalOutputFiles = new HashSet<File>();

    //  private Map<FileDataRecord,DocumentVM> _dataRecordToDataVM = new HashMap<WorkspaceVM.DataRecord, DocumentVM>();

    public DataController(WorkspaceVM<FileDataRecord,?> workspaceViewModel)
    {
        _workspaceViewModel = workspaceViewModel;
    }

    public void addDataFile(final File dataFile) throws IOException
    {

        KnownDatafileFormat format = new DatafileRecogniser().recognise(dataFile);

        final String contents = FileUtils.readFileToString(dataFile);
        final String explorerTitle = dataFile.getName();

        DocumentVM vm = format.execute(new KnownDatafileFormatAlgo<DocumentVM, Void, RuntimeException>() {
            public DocumentVM forSequencings(KnownDatafileFormat format, Void input) throws RuntimeException {
                return new XMLDataVM(contents);
            }

            public DocumentVM forTraces(KnownDatafileFormat format, Void input) throws RuntimeException {
                return new XMLDataVM(contents);
            }

            public DocumentVM forFirstPositive(KnownDatafileFormat format, Void input) throws RuntimeException {
                return new XMLDataVM(contents);
            }

            public DocumentVM forVAALOut(KnownDatafileFormat format, Void input) throws RuntimeException {
                _openVaalOutputFiles.add(dataFile);
                return new VAALOutDataVM(contents);
            }
        }, null);

        final FileDataRecord dataRecord = new FileDataRecord(explorerTitle, dataFile, vm);

        format.execute(new KnownDatafileFormatAlgo<Void, Void, RuntimeException>() {
            public Void forSequencings(KnownDatafileFormat format, Void input) throws RuntimeException {
                _workspaceViewModel.addSequencingsData(dataRecord);
                return null;
            }

            public Void forTraces(KnownDatafileFormat format, Void input) throws RuntimeException {
                _workspaceViewModel.addTraceData(dataRecord);
                return null;
            }

            public Void forFirstPositive(KnownDatafileFormat format, Void input) throws RuntimeException {
                _workspaceViewModel.addFirstPositiveData(dataRecord);
                return null;
            }

            public Void forVAALOut(KnownDatafileFormat format, Void input) throws RuntimeException {
                _workspaceViewModel.addVAALOutData(dataRecord);
                return  null;
            }
        }, null);


        _workspaceViewModel.setFocusDocument(vm);
    }

    public void dataRecordSelected(FileDataRecord dataRecord)
    {
        _workspaceViewModel.setFocusDocument(dataRecord.DocumentVM);
    }

    public void createSequencingDataFileFromOpenVaalOutputFiles() throws IOException, ParseException {
        if(_openVaalOutputFiles.size() < 1)
            return;

        final Map<String,String> sequincingEventIdToVaalOut = new HashMap<String, String>();
        for(File vaalOutFile : _openVaalOutputFiles)
        {
            String id = vaalOutFile.getAbsolutePath();
            if(sequincingEventIdToVaalOut.containsKey(id))
                throw new RuntimeException(); // TODO: do better exception here

            sequincingEventIdToVaalOut.put(id, FileUtils.readFileToString(vaalOutFile));
        }

        File sequencingsXml = new File(_openVaalOutputFiles.iterator().next().getParentFile().getAbsolutePath() + File.separator + "sequencings.xml");
        sequencingsXml.createNewFile();

        String xmlContents = new SequencingsXmlFromVaalOutFactory<String>()
        {
            @Override
            protected String getSourceId(String sequencingEvent) {
                return sequencingEvent;
            }
        }.makeXML(sequincingEventIdToVaalOut);

        FileUtils.writeStringToFile(sequencingsXml, xmlContents);

        addDataFile(sequencingsXml);
    }
}
