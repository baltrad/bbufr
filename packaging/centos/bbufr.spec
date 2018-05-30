%define _prefix /opt/baltrad/bbufr
Name:		bbufr
Version:	2.0.3
Release:	 %{snapshot}%{?dist}
Summary:	BALTRAD interface to EUMETNET OPERA's BUFR software
Group:		Development/Libraries
License:	LGPL
URL:		http://www.baltrad.eu/
Source0:	%{name}-%{version}.tar.gz
BuildRequires:	zlib-devel
BuildRequires:	hdf5-devel
BuildRequires:	proj-devel
BuildRequires:	libpng-devel
Requires:	zlib
Requires:	proj
Requires:	libpng
Requires:	hdf5

%description
BALTRAD interface to EUMETNET OPERA's BUFR software

%prep
%setup -q -n %{name}-%{version}

%build
%{__libtoolize}
%{__aclocal}
%{__autoconf}
%{__automake}
%configure
make

%install
mkdir -p %{buildroot}/opt/baltrad
make install DESTDIR=%{buildroot}

%files
%defattr(-,root,root,-)
%{_prefix}/include/bitio.h
%{_prefix}/include/bufr.h
%{_prefix}/include/bufr_io.h
%{_prefix}/include/bufrlib.h
%{_prefix}/include/desc.h
%{_prefix}/include/rlenc.h
%{_prefix}/%{_lib}/libOperaBufr.a
%{_prefix}/share/bbufr/tables
%changelog
