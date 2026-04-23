from pathlib import Path

from ament_index_python.packages import get_package_share_directory
from launch import LaunchDescription
from launch.actions import (
    DeclareLaunchArgument,
    IncludeLaunchDescription,
    SetEnvironmentVariable,
)
from launch.launch_description_sources import PythonLaunchDescriptionSource
from launch.substitutions import LaunchConfiguration, PathJoinSubstitution


def generate_launch_description() -> LaunchDescription:
    gazebo_share = Path(get_package_share_directory("usv_gazebo"))
    hydro_share = Path(get_package_share_directory("usv_hydro"))
    ros_gz_share = Path(get_package_share_directory("ros_gz_sim"))

    default_world = "calm_sea.sdf"
    default_hydro_config = str(hydro_share / "config" / "hydro_profiles.yaml")

    world_arg = DeclareLaunchArgument(
        "world",
        default_value=default_world,
        description="World file name under usv_gazebo/worlds",
    )
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

    world_path = PathJoinSubstitution(
        [str(gazebo_share / "worlds"), LaunchConfiguration("world")]
    )

    gz = IncludeLaunchDescription(
        PythonLaunchDescriptionSource(str(ros_gz_share / "launch" / "gz_sim.launch.py")),
        launch_arguments={"gz_args": ["-r ", world_path]}.items(),
    )

    return LaunchDescription(
        [
            world_arg,
            hydro_config_arg,
            hydro_profile_arg,
            set_hydro_config_env,
            set_hydro_profile_env,
            gz,
        ]
    )
