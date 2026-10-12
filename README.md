# Robotic-Manipulator-MATLAB-simulation-using-URDF-file
This project simulates the movement in space of a robotic manipulator using the MATLAB program.

# How?
There are no libraries needed to be imported since the goal is to extract the kinematic and dynamic information of the robot through it's URDF file, and then define the forward dynamics equation of the system. 

# More in depth
For defining the system as described above, we calculated analytically the torques, screw axes, rotational axes, Transformation matrices for joint-to-joint cooherence and also the Coriolis matrix and the Mass matrix of the system.

The manipulator consists of five links and five joints where the joint _Arm_link_0_ refers to the base robot. A link is a rigid part of the robot, while a joint connects two links and defines how they can move relative to each other. In other worlds, a joint describes the connection between a parent link and a child link, including its type, origin, and axis when applicable. 

For example, Arm_joint_5 is fixed, meaning the end-effector has no relative motion with respect to the link to which it is attached. Arm_joint_5 is a fixed joint connecting Arm_link_4 to the end-effector. That means the end-effector maintains a fixed position and orientation relative to Arm_link_4.

In this version, we focus more on the practical implementation of the manipulator rather than the full mathematical analysis. Thus, many mathematical stages are skipped in the report. If the reader wishes for the whole project, write in here.
