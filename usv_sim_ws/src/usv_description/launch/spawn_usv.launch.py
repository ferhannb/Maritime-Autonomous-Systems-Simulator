from pathlib import Path

from ament_index_python.packages import get_package_share_directory
from launch import LaunchDescription
from launch.actions import IncludeLaunchDescription, TimerAction
from launch.launch_description_sources import PythonLaunchDescriptionSource
from launch_ros.actions import Node


def generate_launch_description() -> LaunchDescription:
    pkg_share = Path(get_package_share_directory("usv_description"))
    model_path = pkg_share / "sdf" / "usv_minimal.sdf"

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

    return LaunchDescription([gz_launch, TimerAction(period=5.0, actions=[spawn_usv])])
