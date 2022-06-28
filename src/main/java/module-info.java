module com.marginallyclever.makelangelo {
    requires com.marginallyclever.nodegraphcore;
    requires com.formdev.flatlaf;
    requires org.slf4j;
    requires org.json;
    requires org.apache.commons.io;
    requires org.jetbrains.annotations;
    requires org.reflections;
    requires io.sf.carte.echosvg;
    requires java.desktop;
    requires jogamp.fat;
    requires java.prefs;
    requires jssc;
    requires jrpicam;
    requires kabeja;
    requires logback.core;
    requires vecmath;

    opens com.marginallyclever.convenience;
    opens com.marginallyclever.makelangelo.makeart.io.vector;
    opens com.marginallyclever.makelangelo.plotter.plottercontrols;
    opens com.marginallyclever.makelangelo.turtle;

    exports com.marginallyclever.convenience.log to logback.core;

    exports com.marginallyclever.makelangelo.donatelloimpl to com.marginallyclever.nodegraphcore;
    exports com.marginallyclever.makelangelo.donatelloimpl.nodes to com.marginallyclever.nodegraphcore;
    exports com.marginallyclever.makelangelo.donatelloimpl.nodes.shapes to com.marginallyclever.nodegraphcore;

    provides com.marginallyclever.nodegraphcore.NodeRegistry with
            com.marginallyclever.makelangelo.donatelloimpl.DonatelloRegistry;

    provides com.marginallyclever.nodegraphcore.DAORegistry with
            com.marginallyclever.makelangelo.donatelloimpl.DonatelloRegistry;
}