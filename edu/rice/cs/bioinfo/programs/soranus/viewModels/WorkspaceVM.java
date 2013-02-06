package edu.rice.cs.bioinfo.programs.soranus.viewModels;

import edu.rice.cs.bioinfo.library.epidemiology.transmissionMap.Snitkin2012.SnitkinEdge;
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
    class DataRecord
    {
        public final String Title;

        DataRecord(String title) {
            Title = title;
        }
    }

    public final String Title = "Soranus";

    private Set<DataRecord> _sequencingsDataRecords = new HashSet<DataRecord>();

    private Set<DataRecord> _traceDataRecords = new HashSet<DataRecord>();

    private Set<DataRecord> _firstPositiveDataRecords = new HashSet<DataRecord>();

    private Set<Proc1<String>> _dataRecordAddedListeners = new HashSet<Proc1<String>>();

    public void addDataRecordAddedListener(Proc1<String> listener)
    {
        _dataRecordAddedListeners.add(listener);
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


    public void addSequencingsData(String explorerTitle)
    {
        addDataHelp(explorerTitle, _sequencingsDataRecords);
    }

    public void addTraceData(String explorerTitle)
    {
        addDataHelp(explorerTitle, _traceDataRecords);
    }

    public void addFirstPositiveData(String explorerTitle) {
        addDataHelp(explorerTitle, _firstPositiveDataRecords);
    }

    private void addDataHelp(String title, Set<DataRecord> toAdd)
    {
        DataRecord dataRecord = new DataRecord(title);
        toAdd.add(dataRecord);
        for(Proc1<String> dataRecordAddedListener : _dataRecordAddedListeners)
        {
            dataRecordAddedListener.execute(title);
        }
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
