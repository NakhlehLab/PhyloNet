package edu.rice.cs.bioinfo.programs.soranus.models.factories;

import org.w3c.dom.Document;
import org.w3c.dom.Element;
import org.w3c.dom.Node;
import org.w3c.dom.NodeList;
import org.xml.sax.SAXException;

import javax.xml.parsers.DocumentBuilder;
import javax.xml.parsers.DocumentBuilderFactory;
import javax.xml.parsers.ParserConfigurationException;
import java.io.File;
import java.io.IOException;
import java.util.HashMap;
import java.util.Map;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 2/3/13
 * Time: 6:36 PM
 * To change this template use File | Settings | File Templates.
 */
public abstract class PatientTraceMapFactory<I,D,L>
{
    public Map<I, Map<D,L>> makeMapFromXMLFile(File xmlFile) throws ParserConfigurationException, IOException, SAXException {

        DocumentBuilderFactory dbf = DocumentBuilderFactory.newInstance();

         // use the factory to take an instance of the document builder
         DocumentBuilder db = dbf.newDocumentBuilder();
         // parse using the builder to get the DOM mapping of the
         // XML file
         Document dom = db.parse(xmlFile);


         Element doc = dom.getDocumentElement();
         if(doc.getTagName().toLowerCase().equals("traces"))
         {
             Map<I, Map<D,L>> resultMap = new HashMap<I, Map<D,L>>();

             NodeList tracePointElements = dom.getElementsByTagName("TracePoint");
             for(int i = 0; i<tracePointElements.getLength(); i++)
             {
                 Node element = tracePointElements.item(i);
                 String dateText = element.getAttributes().getNamedItem("date").getNodeValue();
                 String entityIdText = element.getAttributes().getNamedItem("entityId").getNodeValue();
                 String locationText = element.getAttributes().getNamedItem("location").getNodeValue();
                 I id = makeId(entityIdText);
                 D date = makeDate(dateText);
                 L location = makeLocation(locationText);

                 if(!resultMap.containsKey(id))
                     resultMap.put(id, new HashMap<D, L>());

                 resultMap.get(id).put(date,location);
             }

             return resultMap;
         }

         throw new IllegalArgumentException();
     }

    protected abstract L makeLocation(String locationText);

    protected abstract D makeDate(String dateText);

    protected abstract I makeId(String sourceIdText) ;

}
