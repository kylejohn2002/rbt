package com.marginallyclever.makelangelo.makeart.tools;

import com.marginallyclever.makelangelo.Translator;
import com.marginallyclever.makelangelo.makeart.TurtleModifierAction;
import com.marginallyclever.makelangelo.turtle.Turtle;
import org.slf4j.Logger;
import org.slf4j.LoggerFactory;

import java.io.Serial;

public class InfillTurtleAction extends TurtleModifierAction {
	private static final Logger logger = LoggerFactory.getLogger(InfillTurtleAction.class);

	@Serial
	private static final long serialVersionUID = -8653065260609614796L;
	
	public InfillTurtleAction() {
		super(Translator.get("InfillTurtleAction.title"));
	}
	
	@Override
	public Turtle run(Turtle t) {
		InfillTurtle infill = new InfillTurtle();
		Turtle result = new Turtle(t);
		try {
			result.add(infill.run(t));
		} catch (Exception ex) {
			logger.error("Failed to infill", ex);
		}
		return result;
	}

}
