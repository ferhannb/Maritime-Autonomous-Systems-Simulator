from pathlib import Path

from ament_index_python.packages import get_package_share_directory
from launch import LaunchDescription
from launch.actions import (
    DeclareLaunchArgument,
    IncludeLaunchDescription,
    SetEnvironmentVariable,
    TimerAction,
)
from launch.launch_description_sources import PythonLaunchDescriptionSource
from launch.substitutions import LaunchConfiguration
from launch_ros.actions import Node


def generate_launch_description() -> LaunchDescription:
    pkg_share = Path(get_package_share_directory("usv_description"))
    hydro_share = Path(get_package_share_directory("usv_hydro"))
    model_path = pkg_share / "sdf" / "usv_minimal.sdf"
    default_hydro_config = str(hydro_share / "config" / "hydro_profiles.yaml")

    hydro_config_arg = DeclareLaunchArgument(
        "hydro_config_file",
        default_value=default_hydro_config,
        description="Path to hydro config YAML",
    )
    hydro_profile_arg = DeclareLaunchArgument(
        "hydro_profile",
        default_value="",
        description=(
            "Optional hydro profile override (e.g. low_damping, stiffness_on). "
            "Empty keeps profile from SDF."
        ),
    )

    set_hydro_config_env = SetEnvironmentVariable(
        name="USV_HYDRO_CONFIG_FILE",
        value=LaunchConfiguration("hydro_config_file"),
    )
    set_hydro_profile_env = SetEnvironmentVariable(
        name="USV_HYDRO_PROFILE",
        value=LaunchConfiguration("hydro_profile"),
    )

    gz_launch = IncludeLaunchDescription(
        PythonLaunchDescriptionSource(
            str(
                Path(get_package_share_directory("ros_gz_sim"))
                / "launch"
                / "gz_sim.launch.py"
            )
        ),
        launch_arguments={"gz_args": "-r empty.sdf"}.items(),
    )

    spawn_usv = Node(
        package="ros_gz_sim",
        executable="create",
        arguments=[
            "-file",
            str(model_path),
            "-name",
            "usv_minimal",
            "-x",
            "0.0",
            "-y",
            "0.0",
            "-z",
            "0.6",
        ],
        output="screen",
    )

    return LaunchDescription(
        [
            hydro_config_arg,
            hydro_profile_arg,
            set_hydro_config_env,
            set_hydro_profile_env,
            gz_launch,
            TimerAction(period=5.0, actions=[spawn_usv]),
        ]
    )
