from conan import ConanFile
from conan.errors import ConanInvalidConfiguration
from conan.tools.build import check_min_cppstd
from conan.tools.files import copy, get
from conan.tools.layout import basic_layout
import os


required_conan_version = ">=2.0"


class PackageConan(ConanFile):
    name = "twofloat"
    description = "C++ library implementing recent double-word (aka double-double) arithmetics"
    license = "MIT"
    url = "https://github.com/conan-io/conan-center-index"
    homepage = "https://github.com/esa-tu-darmstadt/twofloat"
    version = "0.2.0"
    topics = ("double-double", "double-word", "floating-point", "header-only", "arithmetic")
    package_type = "header-library"
    settings = "os", "arch", "compiler", "build_type"
    no_copy_source = True

    def layout(self):
        basic_layout(self, src_folder="src")

    def requirements(self):
        pass

    # same package ID for any package
    def package_id(self):
        self.info.clear()

    def validate(self):
        # Validate the minimum cpp standard supported when installing the package. For C++ projects only
        check_min_cppstd(self, 17)

    def source(self):
        get(self, f"{self.homepage}/archive/refs/tags/v{self.version}.tar.gz", strip_root=True)

    # Suppress warning message about missing build() method when running Conan
    def build(self):
        pass

    # Copy all files to the package folder
    def package(self):
        copy(self, "LICENSE", self.source_folder, os.path.join(self.package_folder, "licenses"))
        # Prefer CMake.install() or similar in case the upstream offers an official method to install the headers.
        copy(self, "*.hpp", os.path.join(self.source_folder, "include"), os.path.join(self.package_folder, "include"))
        