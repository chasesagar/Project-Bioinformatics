from django.urls import path
from .views import HomePageView, CustomPageView
from . import views


urlpatterns = [
    path('',views.HomePageView, name = 'home'),
    path('about/',views.AboutPageView, name = 'about'),
    path('basic-tool/',views.CustomPageView, name = 'basic'),
    path('translate/',views.TranslatePageView, name = 'translate'), 
    path('rev-comp/',views.RevcompPageView, name = 'rev-comp'),
    path('molecular-weight/',views.MolWtPageView, name = 'molwt'),
    path('transcription/',views.TransPageView, name = 'transcribe'),
    path('back-transcription/',views.BackTransPageView, name = 'back-transcribe'),
    path('contact-us/',views.ContactPageView, name = 'contact')
]