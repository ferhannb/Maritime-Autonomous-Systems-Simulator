# USV Simulation Stack  
### A ROS 2-native marine simulation framework for Autonomous Surface Vehicles

[![ROS 2](https://img.shields.io/badge/ROS%202-Jazzy-blue)](#)
[![Gazebo](https://img.shields.io/badge/Gazebo-Harmonic-orange)](#)
[![Platform](https://img.shields.io/badge/Platform-Ubuntu%2024.04-green)](#)
[![Status](https://img.shields.io/badge/Status-Under%20Development-yellow)](#)
[![License](https://img.shields.io/badge/License-MIT-lightgrey)](#)

---

## Table of Contents

- [Overview](#overview)
- [Why This Project Exists](#why-this-project-exists)
- [Manifesto](#manifesto)
- [Vision](#vision)
- [Core Architecture](#core-architecture)
- [System Layers](#system-layers)
- [Development Roadmap](#development-roadmap)
- [Visual Roadmap](#visual-roadmap)
- [Repository Structure](#repository-structure)
- [ROS 2 Interface Philosophy](#ros-2-interface-philosophy)
- [Scenario and Evaluation Strategy](#scenario-and-evaluation-strategy)
- [Validation Philosophy](#validation-philosophy)
- [Engineering Principles](#engineering-principles)
- [Near-Term Milestones](#near-term-milestones)
- [Long-Term Direction](#long-term-direction)
- [What Success Looks Like](#what-success-looks-like)
- [Final Statement](#final-statement)

---

## Overview

This repository aims to build a **modular, physics-aware, ROS 2-integrated simulation ecosystem** for **Autonomous Surface Vehicles (USVs)**.

The objective is not merely to visualize a boat moving in a simulator.  
The objective is to create a simulation framework that supports:

- guidance, navigation, and control development
- thrust allocation and actuator modeling
- sensor simulation and autonomy stack integration
- disturbance-aware testing
- scenario-based evaluation
- system identification feedback loops
- future HIL/SIL and digital twin workflows

This project is designed around a simple idea:

> A useful simulator is not an animation tool.  
> It is an engineering instrument.

---

## Why This Project Exists

Many marine simulation environments fall into one of two categories:

1. **Visually impressive but dynamically shallow**
2. **Academically interesting but difficult to integrate into real robotics workflows**

That gap matters.

For USV development, especially in ROS 2-based autonomy systems, simulation must do more than render a floating vehicle. It must support:

- meaningful control development
- realistic actuator constraints
- reproducible experiments
- measurable performance
- clear separation between truth, estimates, and commands
- eventual comparison with real-world logs

This repository exists to address that need.

---

## Manifesto

We believe a serious USV simulator must be:

- **ROS 2-native**
- **Modular by design**
- **Physics-informed**
- **Control-oriented**
- **Scenario-driven**
- **Validation-ready**
- **Extensible toward real-world deployment workflows**

We do **not** aim to build a monolithic simulator that hides everything behind opaque plugins and fragile dependencies.

We aim to build a **clear, inspectable, replaceable, engineering-first simulation stack**.

---

## Vision

The long-term goal is to establish a simulation framework that can evolve from a research platform into a high-confidence engineering environment.

This means combining:

- **Gazebo Sim** as the simulation world and physics backbone
- **ROS 2** as the communication and autonomy backbone
- **Custom vehicle dynamics and actuator models**
- **Scenario-based evaluation and benchmarking**
- **Sensorized autonomy workflows**
- **Validation against real-world data**
- **A path toward digital twin alignment**

---

## Core Architecture

```text
+------------------------------------------------------------------+
|                    Scenario & Evaluation Layer                   |
| missions | waypoint tasks | DP cases | COLREG cases | metrics    |
+------------------------------------------------------------------+
|                  Navigation, Guidance & Control                  |
| path planner | mission manager | LOS | MPC | allocator | EKF     |
+------------------------------------------------------------------+
|                    Sensor & ROS 2 Interface Layer                |
| IMU | GNSS | heading | AIS | lidar | camera | TF | /clock        |
+------------------------------------------------------------------+
|                 Vehicle Dynamics & Actuator Layer                |
| 3-DOF/6-DOF | thrust model | lag | deadband | saturation         |
+------------------------------------------------------------------+
|                   Marine Physics & World Layer                   |
| buoyancy | hydrodynamics | current | wind | collision            |
+------------------------------------------------------------------+
|                         Gazebo Sim Backbone                      |
+------------------------------------------------------------------+
---
## Visiual Roadmap

[Phase 0]
Foundation
    |
    v
[Phase 1]
Minimal Floating Vehicle
    |
    v
[Phase 2]
Control Integration
    |
    v
[Phase 3]
Actuator Realism
    |
    v
[Phase 4]
Disturbance & Environment
    |
    v
[Phase 5]
Sensorized Autonomy Stack
    |
    v
[Phase 6]
Scenario Library
    |
    v
[Phase 7]
Validation & Digital Twin Direction
