package edu.rice.cs.bioinfo.programs.isolatepipeline;

import org.apache.commons.io.FileUtils;

import java.io.File;
import java.sql.*;
import java.util.List;

/**
 * Created with IntelliJ IDEA.
 * User: Matt
 * Date: 8/4/13
 * Time: 4:33 PM
 * To change this template use File | Settings | File Templates.
 */
public class LoadIsolatesIntoDB
{
    public static void main(String[] args) throws SQLException, ClassNotFoundException, Exception
    {

        Class.forName("com.mysql.jdbc.Driver");

        Connection conn = DriverManager.getConnection("jdbc:mysql://localhost:3306/data", "pipeline", "pipeline");
        PreparedStatement insert = conn.prepareStatement("UPDATE isolates SET collectionDateUtc=? WHERE id=?");

        List<String> lines = FileUtils.readLines(new File("W:\\unorganized\\Logbook_for_Luay.csv"));

        for(String line : lines)
        {
            String[] parts = line.split(",");
            String name = parts[0] + "-" + parts[1];
            String date = parts[2];
            if(date.contains("/"))
            {
                String[] dateparts = date.split("/");
                if(dateparts.length == 3)
                {

                    insert.setString(1, Integer.parseInt(dateparts[2]) + "-" + dateparts[0] + "-" + dateparts[1]);
                    insert.setString(2, name);
                    System.out.println(date);
                    insert.execute();
                }
            }

        }

      /*  PreparedStatement insert = conn.prepareStatement("INSERT INTO isolates (id) VALUES (?)");


        try
        {

            for(File file : new File(args[0]).listFiles())
            {
                if(file.isFile() && file.getName().endsWith(".fastq"))
                {
                    String isolateName = file.getName().substring(0, file.getName().length()-".fastq".length());

                    System.out.println(isolateName);
                    insert.setString(1, isolateName);
                    insert.execute();

                }
            }
        }
        finally {

            insert.close();
        }       */
    }
}
