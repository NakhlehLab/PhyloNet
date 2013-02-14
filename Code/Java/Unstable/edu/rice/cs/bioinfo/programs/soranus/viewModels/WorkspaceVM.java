package edu.rice.cs.bioinfo.programs.soranus.viewModels;

import edu.rice.cs.bioinfo.library.programming.Proc1;

import java.util.HashSet;
import java.util.Set;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 1/30/13
 * Time: 6:11 PM
 * To change this template use File | Settings | File Templates.
 */
public class WorkspaceVM
{
    public class DataRecord
    {
        public final String Title;

        DataRecord(String title) {
            Title = title;
        }
    }

    public class AnalysisRecord
    {
        public final String Title;

        AnalysisRecord(String title)
        {
            Title = title;
        }
    }

    public final String Title = "Soranus";

    private Set<DataRecord> _sequencingsDataRecords = new HashSet<DataRecord>();

    private Set<DataRecord> _traceDataRecords = new HashSet<DataRecord>();

    private Set<DataRecord> _firstPositiveDataRecords = new HashSet<DataRecord>();

    private Set<Proc1<DataRecord>> _dataRecordAddedListeners = new HashSet<Proc1<DataRecord>>();

    public void addDataRecordAddedListener(Proc1<DataRecord> listener)
    {
        _dataRecordAddedListeners.add(listener);
    }

    private Set<Proc1<AnalysisRecord>> _anaylisisAddedListeners = new HashSet<Proc1<AnalysisRecord>>();

    public void addAnalysisAddedListener(Proc1<AnalysisRecord> listener)
    {
        _anaylisisAddedListeners.add(listener);
    }

    private DocumentVM _focusDocument = null;

    public void setFocusDocument(DocumentVM focusDocument)
    {
        _focusDocument = focusDocument;

        for(Proc1<DocumentVM> listener : _focusDocumentChangedListeners)
        {
            listener.execute(_focusDocument);
        }
    }


    private Set<Proc1<DocumentVM>> _focusDocumentChangedListeners = new HashSet<Proc1<DocumentVM>>();

    public void addFocusDocumentChangedListener(Proc1<DocumentVM> listener)
    {
        _focusDocumentChangedListeners.add(listener);
    }


    public DataRecord addSequencingsData(String explorerTitle)
    {
        return addDataHelp(explorerTitle, _sequencingsDataRecords);
    }

    public DataRecord addTraceData(String explorerTitle)
    {
        return addDataHelp(explorerTitle, _traceDataRecords);
    }

    public DataRecord addFirstPositiveData(String explorerTitle) {
        return addDataHelp(explorerTitle, _firstPositiveDataRecords);
    }

    private DataRecord addDataHelp(String title, Set<DataRecord> toAdd)
    {
        DataRecord dataRecord = new DataRecord(title);
        toAdd.add(dataRecord);
        for(Proc1<DataRecord> dataRecordAddedListener : _dataRecordAddedListeners)
        {
            dataRecordAddedListener.execute(dataRecord);
        }

        return dataRecord;
    }

    public AnalysisRecord addAnalysis(String title)
    {
        AnalysisRecord record = new AnalysisRecord(title);
        for(Proc1<AnalysisRecord> listener : _anaylisisAddedListeners)
        {
            listener.execute(record);
        }

        return  record;
    }

    public boolean isDataRecordTitleRepresentingSequencingsData(String title)
    {
        for(DataRecord dr : _sequencingsDataRecords)
        {
            if(dr.Title.equals(title))
            {
                return true;
            }
        }

        return false;
    }

    public boolean isDataRecordTitleRepresentingTraceData(String title)
    {
        for(DataRecord dr : _traceDataRecords)
        {
            if(dr.Title.equals(title))
            {
                return true;
            }
        }

        return false;
    }

    public boolean isDataRecordTitleRepresentingFirstPositiveData(String title)
    {
        for(DataRecord dr : _firstPositiveDataRecords)
        {
            if(dr.Title.equals(title))
            {
                return true;
            }
        }

        return false;
    }

}
